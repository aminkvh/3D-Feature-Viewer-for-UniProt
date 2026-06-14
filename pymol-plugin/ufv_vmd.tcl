# =============================================================================================
# 3D Feature Viewer for UniProt - VMD plugin
# =============================================================================================
# Projects UniProt residue-level annotations (PTMs, disease variants, ClinVar/AlphaMissense,
# functional sites, domains/regions, membrane topology) onto a structure or trajectory loaded in
# VMD - or downloads the AlphaFold model directly.
#
# Fetching/mapping is delegated to the bundled Python backend (ufv_pymol.py), which uses the same
# public data sources as the browser extension (UniProtKB, EBI Proteins API, AlphaFold DB,
# PDBe/SIFTS). All projection is done natively in VMD via atomselect + representations.
#
# Requirements: python3 on PATH (set another with `ufv_python <path>`), and ufv_pymol.py next to
# this file (or set with `ufv_backend <path>`).
#
# Install:  source ufv_vmd.tcl        (or place in a VMD plugin dir)
#
# Quick start
#   ufv_load P35498                 ;# download AlphaFold model + annotate
#   ufv_gui                         ;# graphical panel
#
# Annotate a structure/trajectory you already loaded (molid defaults to top):
#   mol new mytraj.pdb
#   ufv_fetch P35498
#   ufv_map identity                ;# resid == UniProt position
#   ufv_map sifts 7dtd              ;# map through PDBe/SIFTS for PDB 7DTD
#   ufv_chain A 200 5 480          ;# chain A resid 5 == UniProt 200, valid 5..480
#   ufv_ptms ; ufv_variants ; ufv_sites ; ufv_domains ; ufv_topology ; ufv_alphamissense
#   ufv_clear
# =============================================================================================

package provide ufv 1.0

namespace eval ::ufv {
    variable python  "python"
    variable backend [file join [file dirname [file normalize [info script]]] ufv_pymol.py]
    variable uid     ""
    variable molid   ""
    variable mode    "identity"
    variable colorid 33
    variable manual
    array set manual {}
    # annotation arrays a_* / sifts s_* are created by sourcing the backend output
}

# ---- python / backend configuration --------------------------------------------------------
proc ufv_python  {path} { set ::ufv::python  $path; puts "\[UFV\] python = $path" }
proc ufv_backend {path} { set ::ufv::backend $path; puts "\[UFV\] backend = $path" }

proc ::ufv::run {args} {
    variable python
    variable backend
    if {![file exists $backend]} {
        error "backend not found: $backend  (set with: ufv_backend <path to ufv_pymol.py>)"
    }
    return [exec $python $backend {*}$args]
}

# ---- fetch / map ----------------------------------------------------------------------------
proc ::ufv::fetch {uid} {
    variable uid
    # Clear any previous annotation arrays so stale layers don't linger.
    foreach v [info vars ::ufv::a_*] { unset -nocomplain $v }
    eval [::ufv::run --emit-tcl $uid]
    set ::ufv::uid $uid
    puts "[format {[UFV] %s: %d PTMs, %d variants, %d sites, %d domains, %d topology} \
        $uid [llength $::ufv::a_ptms] [llength $::ufv::a_variants] [llength $::ufv::a_sites] \
        [llength $::ufv::a_domains] [llength $::ufv::a_topology]]"
}

proc ::ufv::need {uid} {
    if {$uid ne ""} { if {$uid ne $::ufv::uid} { ::ufv::fetch $uid } ; return }
    if {$::ufv::uid eq ""} { error "no accession fetched yet - run: ufv_fetch <uniprot_id>" }
}

proc ufv_fetch {uid} { ::ufv::fetch $uid }

proc ufv_map {args} {
    # ufv_map identity            |  ufv_map sifts <pdbid>
    set m [lindex $args 0]
    if {$m eq "sifts"} {
        set pdb [lindex $args 1]
        if {$::ufv::uid eq ""} { error "fetch an accession first: ufv_fetch <uniprot_id>" }
        foreach v [info vars ::ufv::s_*] { unset -nocomplain $v }
        eval [::ufv::run --emit-sifts $pdb $::ufv::uid]
        set ::ufv::mode sifts
        puts "\[UFV\] numbering: SIFTS via [string toupper $pdb] (chains $::ufv::s_chains)"
    } else {
        set ::ufv::mode identity
        puts "\[UFV\] numbering: identity (resid == UniProt position)"
    }
}

proc ufv_chain {chain ustart {rstart {}} {rend {}}} {
    if {$rstart eq ""} { set rstart $ustart }
    set ::ufv::manual($chain) [list $ustart $rstart $rend]
    set ::ufv::mode manual
    puts "\[UFV\] chain $chain: resid $rstart == UniProt $ustart[expr {$rend ne {} ? \" (..$rend)\" : {}}]"
}

# ---- numbering ------------------------------------------------------------------------------
proc ::ufv::uni_to_resi {chain pos} {
    switch $::ufv::mode {
        identity { return $pos }
        manual {
            variable manual
            if {![info exists manual($chain)]} { return "" }
            lassign $manual($chain) us rs re
            set resi [expr {$rs + ($pos - $us)}]
            if {$resi < $rs} { return "" }
            if {$re ne "" && $resi > $re} { return "" }
            return $resi
        }
        sifts {
            set var ::ufv::s_$chain
            if {![info exists $var]} { return "" }
            foreach seg [set $var] {
                lassign $seg us ue ps pe
                if {$pos >= $us && $pos <= $ue} { return [expr {$ps + ($pos - $us)}] }
            }
            return ""
        }
    }
    return $pos
}

proc ::ufv::mapped_chains {molid} {
    switch $::ufv::mode {
        sifts  { return $::ufv::s_chains }
        manual { return [array names ::ufv::manual] }
        default {
            set s [atomselect $molid "protein"]
            set c [lsort -unique [$s get chain]]
            $s delete
            if {[llength $c] == 0} { return [list ""] }
            return $c
        }
    }
}

# ---- helpers --------------------------------------------------------------------------------
proc ::ufv::resolve {molid} {
    if {$molid eq "" || $molid eq "top"} {
        if {$::ufv::molid ne "" && [lsearch [molinfo list] $::ufv::molid] >= 0} { return $::ufv::molid }
        return [molinfo top get id]
    }
    return $molid
}

proc ::ufv::next_color {hex} {
    variable colorid
    set h [string trimleft $hex "#"]
    scan [string range $h 0 1] %x r
    scan [string range $h 2 3] %x g
    scan [string range $h 4 5] %x b
    set id $colorid
    color change rgb $id [expr {$r/255.0}] [expr {$g/255.0}] [expr {$b/255.0}]
    incr colorid
    if {$colorid > 1000} { set colorid 33 }
    return $id
}

# Build a VMD selection string for a list of UniProt positions on the resolved molecule.
proc ::ufv::sel_positions {molid positions caonly} {
    set clauses {}
    foreach ch [::ufv::mapped_chains $molid] {
        set resids {}
        foreach p $positions {
            set r [::ufv::uni_to_resi $ch $p]
            if {$r ne ""} { lappend resids $r }
        }
        if {[llength $resids]} {
            set resids [lsort -integer -unique $resids]
            if {$ch eq ""} {
                lappend clauses "(resid [join $resids { }])"
            } else {
                lappend clauses "(chain $ch and resid [join $resids { }])"
            }
        }
    }
    if {![llength $clauses]} { return "" }
    set sel "([join $clauses { or }])"
    if {$caonly} { set sel "$sel and name CA" }
    return $sel
}

proc ::ufv::addrep {molid style sel colorid} {
    mol representation $style
    mol selection $sel
    if {$colorid ne ""} { mol color ColorID $colorid } else { mol color Name }
    mol material Opaque
    mol addrep $molid
}

# ---- projection layers ----------------------------------------------------------------------
proc ::ufv::point_layer {molid items tag} {
    # items: list of {pos color [token]}; draws coloured VDW spheres on Ca, grouped by colour.
    array set groups {}
    foreach e $items {
        set pos   [lindex $e 0]
        set color [lindex $e 1]
        lappend groups($color) $pos
    }
    set n 0
    foreach color [array names groups] {
        set sel [::ufv::sel_positions $molid $groups($color) 1]
        if {$sel eq ""} continue
        ::ufv::addrep $molid "VDW 0.5 12" $sel [::ufv::next_color $color]
        incr n
    }
    array unset groups
    return $n
}

proc ufv_ptms {{molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    set n [::ufv::point_layer $molid $::ufv::a_ptms ptm]
    puts "\[UFV\] PTMs projected ($n colour groups)."
}

proc ufv_variants {args} {
    # ufv_variants [molid] [filter]   filter: pathogenic|benign|uncertain|deleterious
    set molid ""; set filter ""
    foreach a $args {
        if {[lsearch {pathogenic benign uncertain deleterious} $a] >= 0} { set filter $a } else { set molid $a }
    }
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    set items {}
    foreach e $::ufv::a_variants {
        if {$filter ne "" && [lindex $e 2] ne $filter} continue
        lappend items $e
    }
    set n [::ufv::point_layer $molid $items var]
    puts "\[UFV\] variants projected ($n colour groups[expr {$filter ne {} ? \", $filter only\" : {}}])."
}

proc ufv_sites {{molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    set n [::ufv::point_layer $molid $::ufv::a_sites site]
    puts "\[UFV\] functional sites projected ($n)."
}

proc ufv_domains {{molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    set n 0
    foreach e $::ufv::a_domains {
        lassign $e start end color isrange
        set positions {}
        for {set p $start} {$p <= $end} {incr p} { lappend positions $p }
        if {$isrange} {
            set sel [::ufv::sel_positions $molid $positions 0]
            if {$sel eq ""} continue
            ::ufv::addrep $molid "NewCartoon" $sel [::ufv::next_color $color]
        } else {
            set sel [::ufv::sel_positions $molid $positions 1]
            if {$sel eq ""} continue
            ::ufv::addrep $molid "VDW 0.6 12" $sel [::ufv::next_color $color]
        }
        incr n
    }
    puts "\[UFV\] $n domain/region features projected."
}

proc ufv_topology {{molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    set n 0
    foreach e $::ufv::a_topology {
        lassign $e start end color
        set positions {}
        for {set p $start} {$p <= $end} {incr p} { lappend positions $p }
        set sel [::ufv::sel_positions $molid $positions 0]
        if {$sel eq ""} continue
        ::ufv::addrep $molid "NewCartoon" $sel [::ufv::next_color $color]
        incr n
    }
    puts "\[UFV\] $n topology segments projected."
}

proc ufv_alphamissense {{molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    if {[llength $::ufv::a_am] == 0} { puts "\[UFV\] no AlphaMissense scores."; return }
    array set buckets {}
    foreach e $::ufv::a_am {
        lassign $e pos avg
        if {$avg >= 0.78}      { set c "#b71c1c" } \
        elseif {$avg >= 0.564} { set c "#e06666" } \
        elseif {$avg >= 0.34}  { set c "#b9c2cf" } \
        else                   { set c "#3d85c8" }
        lappend buckets($c) $pos
    }
    foreach c [array names buckets] {
        set sel [::ufv::sel_positions $molid $buckets($c) 0]
        if {$sel ne ""} { ::ufv::addrep $molid "NewCartoon" $sel [::ufv::next_color $c] }
    }
    array unset buckets
    puts "\[UFV\] AlphaMissense colouring applied."
}

proc ufv_clear {{molid ""}} {
    set molid [::ufv::resolve $molid]
    set n [molinfo $molid get numreps]
    for {set i [expr {$n - 1}]} {$i >= 0} {incr i -1} { mol delrep $i $molid }
    mol representation NewCartoon
    mol selection {protein}
    mol color ColorID 6
    mol material Opaque
    mol addrep $molid
    puts "\[UFV\] cleared overlays on mol $molid."
}

# ---- load (download AlphaFold model) --------------------------------------------------------
proc ufv_load {uid} {
    set path [string trim [::ufv::run --download $uid]]
    if {$path eq "" || ![file exists $path]} { puts "\[UFV\] download failed for $uid"; return }
    set m [mol new $path waitfor all]
    set ::ufv::molid $m
    mol delrep 0 $m
    mol representation NewCartoon
    mol selection {protein}
    mol color ColorID 6
    mol material Opaque
    mol addrep $m
    ::ufv::fetch $uid
    set ::ufv::mode identity
    ufv_ptms $m
    ufv_variants $m
    ufv_sites $m
    puts "\[UFV\] loaded and annotated $uid as mol $m. Try: ufv_domains | ufv_topology | ufv_alphamissense"
}

# ---- minimal Tk GUI -------------------------------------------------------------------------
proc ufv_gui {} {
    if {[winfo exists .ufv]} { wm deiconify .ufv; raise .ufv; return }
    toplevel .ufv
    wm title .ufv "3D Feature Viewer for UniProt"

    set f [frame .ufv.f -padx 10 -pady 10]; pack $f -fill both -expand 1

    label $f.lu -text "UniProt accession:"; grid $f.lu -row 0 -column 0 -sticky w
    entry $f.eu -width 14 -textvariable ::ufv::gui_uid; grid $f.eu -row 0 -column 1 -sticky w
    set ::ufv::gui_uid $::ufv::uid

    labelframe $f.num -text "Residue numbering" -padx 6 -pady 6
    grid $f.num -row 1 -column 0 -columnspan 2 -sticky ew -pady 6
    radiobutton $f.num.i -text "Identity (resid == UniProt)" -variable ::ufv::gui_mode -value identity
    grid $f.num.i -row 0 -column 0 -columnspan 3 -sticky w
    radiobutton $f.num.s -text "SIFTS, PDB id:" -variable ::ufv::gui_mode -value sifts
    entry $f.num.sp -width 8 -textvariable ::ufv::gui_pdb
    grid $f.num.s -row 1 -column 0 -sticky w; grid $f.num.sp -row 1 -column 1 -sticky w
    radiobutton $f.num.m -text "Manual chain:" -variable ::ufv::gui_mode -value manual
    grid $f.num.m -row 2 -column 0 -sticky w
    frame $f.num.mf; grid $f.num.mf -row 2 -column 1 -columnspan 2 -sticky w
    foreach {lbl var w} {chain ::ufv::gui_ch 3 U@ ::ufv::gui_us 6 resid@ ::ufv::gui_rs 5 end ::ufv::gui_re 5} {
        label $f.num.mf.l$var -text $lbl; entry $f.num.mf.e$var -width $w -textvariable $var
        pack $f.num.mf.l$var $f.num.mf.e$var -side left
    }
    set ::ufv::gui_mode identity; set ::ufv::gui_ch A
    button $f.num.apply -text "Apply numbering" -command {
        if {$::ufv::gui_uid eq ""} { tk_messageBox -message "Enter an accession."; return }
        if {$::ufv::uid ne $::ufv::gui_uid} { ufv_fetch $::ufv::gui_uid }
        switch $::ufv::gui_mode {
            sifts  { ufv_map sifts $::ufv::gui_pdb }
            manual { ufv_chain $::ufv::gui_ch $::ufv::gui_us $::ufv::gui_rs $::ufv::gui_re }
            default { ufv_map identity }
        }
    }
    grid $f.num.apply -row 3 -column 0 -columnspan 3 -sticky w -pady 4

    labelframe $f.lay -text "Annotation layers" -padx 6 -pady 6
    grid $f.lay -row 2 -column 0 -columnspan 2 -sticky ew
    set i 0
    foreach {lbl cmd} {PTMs ufv_ptms "Disease variants" ufv_variants "Functional sites" ufv_sites \
                       "Domains / regions" ufv_domains "Membrane topology" ufv_topology \
                       AlphaMissense ufv_alphamissense} {
        button $f.lay.b$i -text $lbl -width 18 -command "if {\$::ufv::gui_uid eq {}} { tk_messageBox -message {Enter an accession.} } else { ufv_fetch \$::ufv::gui_uid ; $cmd }"
        grid $f.lay.b$i -row [expr {$i/2}] -column [expr {$i%2}] -sticky w -padx 2 -pady 2
        incr i
    }

    frame $f.bot; grid $f.bot -row 3 -column 0 -columnspan 2 -sticky w -pady 6
    button $f.bot.load -text "Download AlphaFold + annotate" -command { ufv_load $::ufv::gui_uid }
    button $f.bot.clear -text "Clear" -command { ufv_clear }
    pack $f.bot.load $f.bot.clear -side left -padx 2
}

# VMD menu registration (when loaded as an extension).
if {![catch {package present vmd}]} {
    catch { vmd_install_extension ufv ufv_gui "Analysis/3D Feature Viewer for UniProt" }
}

puts "\[UFV\] 3D Feature Viewer for UniProt (VMD) loaded. Commands: ufv_load, ufv_fetch, ufv_map,"
puts "      ufv_chain, ufv_ptms, ufv_variants, ufv_sites, ufv_domains, ufv_topology,"
puts "      ufv_alphamissense, ufv_clear, ufv_gui  (set interpreter with: ufv_python <path>)"
