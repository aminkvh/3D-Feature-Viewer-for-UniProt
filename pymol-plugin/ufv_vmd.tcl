# =============================================================================================
# 3D Feature Viewer for UniProt - VMD plugin
# =============================================================================================
# Projects UniProt residue-level annotations (PTMs, disease variants, ClinVar/AlphaMissense,
# functional sites, domains/regions, membrane topology, mutagenesis, mutation burden) onto a
# structure or trajectory loaded in VMD - or downloads/loads AlphaFold / experimental / computed
# models directly. A residue readout shows the per-residue report (variants + disease + gnomAD and
# ProtVar predictors EVE/ESM1b/conservation/FoldX) on demand.
#
# Fetching/mapping/prediction is delegated to the bundled Python backend (ufv_pymol.py), which uses
# the same public data sources as the browser extension and PyMOL plugin (UniProtKB, EBI Proteins
# API, AlphaFold DB, PDBe/SIFTS, EMBL-EBI ProtVar). All projection is done natively in VMD via
# atomselect + representations.
#
# Requirements: python3 on PATH (set another with `ufv_python <path>`), and ufv_pymol.py next to
# this file (or set with `ufv_backend <path>`).
#
# Install:  source ufv_vmd.tcl        (or place in a VMD plugin dir);  then:  ufv_gui
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
    # annotation lists a_* / sifts s_* are created by sourcing the backend output
    variable a_structures {}
    variable a_function   ""
    variable a_protnlm    ""
    variable a_residue    ""
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

proc ::ufv::status {msg} {
    if {[winfo exists .ufv.status]} { .ufv.status configure -text $msg }
    puts "\[UFV\] $msg"
}

# ---- fetch / map ----------------------------------------------------------------------------
proc ::ufv::fetch {uid} {
    variable uid
    foreach v [info vars ::ufv::a_*] { unset -nocomplain $v }
    eval [::ufv::run --emit-tcl $uid]
    set ::ufv::uid $uid
    puts "[format {[UFV] %s: %d PTMs, %d variants, %d sites, %d domains, %d topology, %d mutagenesis} \
        $uid [llength $::ufv::a_ptms] [llength $::ufv::a_variants] [llength $::ufv::a_sites] \
        [llength $::ufv::a_domains] [llength $::ufv::a_topology] [llength $::ufv::a_mutagenesis]]"
}

proc ::ufv::need {uid} {
    if {$uid ne ""} { if {$uid ne $::ufv::uid} { ::ufv::fetch $uid } ; return }
    if {$::ufv::uid eq ""} { error "no accession fetched yet - run: ufv_fetch <uniprot_id>" }
}

proc ufv_fetch {uid} { ::ufv::fetch $uid }

proc ufv_map {args} {
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

proc ::ufv::addrep {molid style sel colorid {material Opaque}} {
    mol representation $style
    mol selection $sel
    if {$colorid ne ""} { mol color ColorID $colorid } else { mol color Name }
    mol material $material
    mol addrep $molid
}

# ---- projection layers ----------------------------------------------------------------------
proc ::ufv::point_layer {molid items tag} {
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
    ::ufv::status "PTMs projected ([::ufv::point_layer $molid $::ufv::a_ptms ptm] groups)."
}

proc ufv_variants {args} {
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
    ::ufv::status "variants projected ([::ufv::point_layer $molid $items var] groups[expr {$filter ne {} ? \", $filter\" : {}}])."
}

proc ufv_sites {{molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    ::ufv::status "functional sites projected ([::ufv::point_layer $molid $::ufv::a_sites site])."
}

proc ufv_mutagenesis {{molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    ::ufv::status "mutagenesis projected ([::ufv::point_layer $molid $::ufv::a_mutagenesis mutag])."
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
    ::ufv::status "$n domain/region features projected."
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
    ::ufv::status "$n topology segments projected."
}

proc ufv_alphamissense {{molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    if {[llength $::ufv::a_am] == 0} { ::ufv::status "no AlphaMissense scores."; return }
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
    ::ufv::status "AlphaMissense colouring applied."
}

proc ufv_burden {{molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    if {[llength $::ufv::a_burden] == 0} { ::ufv::status "no burden-positive residues."; return }
    set sel [::ufv::sel_positions $molid $::ufv::a_burden 0]
    if {$sel ne ""} { ::ufv::addrep $molid "NewCartoon" $sel [::ufv::next_color "#e65100"] }
    ::ufv::status "mutation burden projected ([llength $::ufv::a_burden] residues)."
}

# pLDDT (predicted) / B-factor (experimental): colour the cartoon by the Beta column (0-100).
proc ::ufv::beta_color {molid label} {
    set molid [::ufv::resolve $molid]
    mol representation NewCartoon
    mol selection {protein}
    mol color Beta
    mol material Opaque
    mol addrep $molid
    set rep [expr {[molinfo $molid get numreps] - 1}]
    catch { mol scaleminmax $molid $rep 0 100 }
    ::ufv::status "coloured by Beta column ($label)."
}
proc ufv_plddt   {{molid ""}} { ::ufv::beta_color $molid "pLDDT - predicted models" }
proc ufv_bfactor {{molid ""}} { ::ufv::beta_color $molid "B-factor - experimental" }

proc ufv_clear {{molid ""}} {
    set molid [::ufv::resolve $molid]
    set n [molinfo $molid get numreps]
    for {set i [expr {$n - 1}]} {$i >= 0} {incr i -1} { mol delrep $i $molid }
    mol representation NewCartoon
    mol selection {protein}
    mol color ColorID 6
    mol material Opaque
    mol addrep $molid
    ::ufv::status "cleared overlays on mol $molid."
}

# ---- residue report -------------------------------------------------------------------------
proc ufv_residue {pos {molid ""}} {
    ::ufv::need ""
    catch { unset ::ufv::a_residue }
    eval [::ufv::run --residue $::ufv::uid $pos]
    set txt [expr {[info exists ::ufv::a_residue] ? $::ufv::a_residue : "no data"}]
    # Highlight the residue in 3D (licorice) so it's locatable.
    set molid [::ufv::resolve $molid]
    set sel [::ufv::sel_positions $molid [list $pos] 0]
    if {$sel ne ""} { ::ufv::addrep $molid "Licorice 0.3 12" "byres ($sel)" [::ufv::next_color "#ffa726"] }
    if {[winfo exists .ufv.rep.t]} {
        .ufv.rep.t configure -state normal
        .ufv.rep.t delete 1.0 end
        .ufv.rep.t insert end $txt
        .ufv.rep.t configure -state disabled
    } else {
        puts $txt
    }
    return $txt
}

# ---- structures: list + load ----------------------------------------------------------------
proc ::ufv::fetch_structures {uid} {
    set ::ufv::a_structures {}
    eval [::ufv::run --structures $uid]
    return $::ufv::a_structures
}

proc ::ufv::base_cartoon {m} {
    mol delrep 0 $m
    mol representation NewCartoon
    mol selection {protein}
    mol color ColorID 6
    mol material Opaque
    mol addrep $m
}

proc ufv_load_structure {key} {
    foreach s $::ufv::a_structures {
        if {[lindex $s 0] ne $key} continue
        set url [lindex $s 3]; set fmt [lindex $s 4]; set numbering [lindex $s 5]; set pdb [lindex $s 6]
        ::ufv::status "downloading $key ..."
        set path [string trim [::ufv::run --get $url $fmt]]
        if {$path eq "" || ![file exists $path]} { ::ufv::status "download failed for $key"; return }
        set m [mol new $path waitfor all]
        set ::ufv::molid $m
        ::ufv::base_cartoon $m
        if {$numbering eq "sifts" && $pdb ne "-"} { ufv_map sifts $pdb } else { ufv_map identity }
        ::ufv::status "loaded $key as mol $m."
        return $m
    }
    ::ufv::status "structure $key not in the list - fetch first."
}

# ---- load (download AlphaFold model directly) -----------------------------------------------
proc ufv_load {uid} {
    set path [string trim [::ufv::run --download $uid]]
    if {$path eq "" || ![file exists $path]} { ::ufv::status "download failed for $uid"; return }
    set m [mol new $path waitfor all]
    set ::ufv::molid $m
    ::ufv::base_cartoon $m
    ::ufv::fetch $uid
    set ::ufv::mode identity
    ufv_ptms $m ; ufv_variants $m ; ufv_sites $m
    ::ufv::status "loaded and annotated $uid as mol $m."
}

# ---- Tk GUI ---------------------------------------------------------------------------------
proc ufv_gui {} {
    if {[winfo exists .ufv]} { wm deiconify .ufv; raise .ufv; return }
    toplevel .ufv
    wm title .ufv "3D Feature Viewer for UniProt"

    label .ufv.hdr -text "3D Feature Viewer for UniProt" -font {Helvetica 12 bold} \
        -fg white -bg "#00695c" -anchor w -padx 10 -pady 6
    pack .ufv.hdr -fill x

    set f [frame .ufv.f -padx 10 -pady 8]; pack $f -fill both -expand 1

    # --- accession + fetch ---
    set s1 [labelframe $f.s1 -text "Structure" -padx 6 -pady 6]; pack $s1 -fill x
    label $s1.lu -text "UniProt"; grid $s1.lu -row 0 -column 0 -sticky w
    entry $s1.eu -width 12 -textvariable ::ufv::gui_uid; grid $s1.eu -row 0 -column 1 -sticky w
    set ::ufv::gui_uid $::ufv::uid
    button $s1.fetch -text "Fetch" -command {
        if {$::ufv::gui_uid eq ""} { tk_messageBox -message "Enter an accession."; return }
        ::ufv::status "fetching $::ufv::gui_uid ..."
        ufv_fetch $::ufv::gui_uid
        ::ufv::fetch_structures $::ufv::gui_uid
        .ufv.f.s1.lb delete 0 end
        foreach st $::ufv::a_structures { .ufv.f.s1.lb insert end "[lindex $st 1]  ([lindex $st 2])" }
        ::ufv::ui_function
        ::ufv::status "fetched $::ufv::gui_uid - choose a structure and Load."
    }
    grid $s1.fetch -row 0 -column 2 -sticky w -padx 4

    listbox $s1.lb -height 5 -width 46 -yscrollcommand "$s1.sb set"
    scrollbar $s1.sb -command "$s1.lb yview"
    grid $s1.lb -row 1 -column 0 -columnspan 3 -sticky ew -pady 4
    grid $s1.sb -row 1 -column 3 -sticky ns
    button $s1.load -text "Load selected" -command {
        set i [.ufv.f.s1.lb curselection]
        if {$i eq ""} { tk_messageBox -message "Select a structure from the list."; return }
        ufv_load_structure [lindex [lindex $::ufv::a_structures $i] 0]
    }
    button $s1.loadaf -text "Quick AlphaFold" -command { if {$::ufv::gui_uid ne ""} { ufv_load $::ufv::gui_uid } }
    grid $s1.load -row 2 -column 0 -columnspan 2 -sticky w -pady 2
    grid $s1.loadaf -row 2 -column 2 -columnspan 2 -sticky w -pady 2

    # --- numbering ---
    set num [labelframe $f.num -text "Residue numbering" -padx 6 -pady 6]; pack $num -fill x -pady 6
    radiobutton $num.i -text "Identity (resid == UniProt)" -variable ::ufv::gui_mode -value identity
    grid $num.i -row 0 -column 0 -columnspan 3 -sticky w
    radiobutton $num.s -text "SIFTS, PDB id:" -variable ::ufv::gui_mode -value sifts
    entry $num.sp -width 8 -textvariable ::ufv::gui_pdb
    grid $num.s -row 1 -column 0 -sticky w; grid $num.sp -row 1 -column 1 -sticky w
    radiobutton $num.m -text "Manual chain:" -variable ::ufv::gui_mode -value manual
    grid $num.m -row 2 -column 0 -sticky w
    frame $num.mf; grid $num.mf -row 2 -column 1 -columnspan 2 -sticky w
    foreach {lbl var w} {chain ::ufv::gui_ch 3 U@ ::ufv::gui_us 6 resid@ ::ufv::gui_rs 5 end ::ufv::gui_re 5} {
        label $num.mf.l$var -text $lbl; entry $num.mf.e$var -width $w -textvariable $var
        pack $num.mf.l$var $num.mf.e$var -side left
    }
    set ::ufv::gui_mode identity; set ::ufv::gui_ch A
    button $num.apply -text "Apply numbering" -command {
        if {$::ufv::gui_uid eq ""} { tk_messageBox -message "Enter an accession."; return }
        if {$::ufv::uid ne $::ufv::gui_uid} { ufv_fetch $::ufv::gui_uid }
        switch $::ufv::gui_mode {
            sifts  { ufv_map sifts $::ufv::gui_pdb }
            manual { ufv_chain $::ufv::gui_ch $::ufv::gui_us $::ufv::gui_rs $::ufv::gui_re }
            default { ufv_map identity }
        }
    }
    grid $num.apply -row 3 -column 0 -columnspan 3 -sticky w -pady 4

    # --- marker layers ---
    set lay [labelframe $f.lay -text "Layers (markers)" -padx 6 -pady 6]; pack $lay -fill x
    set i 0
    foreach {lbl cmd} {PTMs ufv_ptms "Disease variants" ufv_variants "Functional sites" ufv_sites \
                       "Mutagenesis" ufv_mutagenesis} {
        button $lay.b$i -text $lbl -width 17 \
            -command "if {\$::ufv::gui_uid eq {}} { tk_messageBox -message {Enter an accession.} } else { ufv_fetch \$::ufv::gui_uid ; $cmd }"
        grid $lay.b$i -row [expr {$i/2}] -column [expr {$i%2}] -sticky w -padx 2 -pady 2
        incr i
    }

    # --- cartoon colouring ---
    set col [labelframe $f.col -text "Cartoon colouring" -padx 6 -pady 6]; pack $col -fill x -pady 6
    set i 0
    foreach {lbl cmd} {Domains ufv_domains Topology ufv_topology AlphaMissense ufv_alphamissense \
                       Burden ufv_burden "pLDDT (predicted)" ufv_plddt "B-factor (exp.)" ufv_bfactor} {
        button $col.b$i -text $lbl -width 17 \
            -command "if {\$::ufv::gui_uid eq {}} { tk_messageBox -message {Enter an accession.} } else { ufv_fetch \$::ufv::gui_uid ; $cmd }"
        grid $col.b$i -row [expr {$i/2}] -column [expr {$i%2}] -sticky w -padx 2 -pady 2
        incr i
    }

    # --- residue report ---
    set rep [labelframe $f.rep -text "Residue report" -padx 6 -pady 6]; pack $rep -fill both -expand 1 -pady 6
    label $rep.l -text "UniProt position:"; grid $rep.l -row 0 -column 0 -sticky w
    entry $rep.e -width 8 -textvariable ::ufv::gui_pos; grid $rep.e -row 0 -column 1 -sticky w
    button $rep.show -text "Show" -command {
        if {$::ufv::gui_pos eq ""} { tk_messageBox -message "Enter a residue position."; return }
        ::ufv::status "fetching residue $::ufv::gui_pos (incl. ProtVar) ..."
        ufv_residue $::ufv::gui_pos
        ::ufv::status "residue $::ufv::gui_pos report shown."
    }
    grid $rep.show -row 0 -column 2 -sticky w -padx 4
    text $rep.t -height 9 -width 52 -wrap word -state disabled -font {Courier 9}
    grid $rep.t -row 1 -column 0 -columnspan 3 -sticky nsew -pady 4
    grid rowconfigure $rep 1 -weight 1
    grid columnconfigure $rep 2 -weight 1

    # --- function context + clear ---
    label $f.func -textvariable ::ufv::a_function -anchor w -justify left -wraplength 380 -fg "#5b3a8e"
    pack $f.func -fill x
    label .ufv.status -text "Ready." -anchor w -relief sunken -padx 6
    pack .ufv.status -fill x -side bottom
    frame $f.bot; pack $f.bot -fill x -pady 4
    button $f.bot.clear -text "Clear overlays" -command { ufv_clear }
    pack $f.bot.clear -side left -padx 2
}

# Refresh the function-context label after a fetch.
proc ::ufv::ui_function {} {
    set t ""
    if {[info exists ::ufv::a_protnlm] && $::ufv::a_protnlm ne ""} { append t $::ufv::a_protnlm "\n" }
    if {[info exists ::ufv::a_function] && $::ufv::a_function ne ""} { append t $::ufv::a_function }
    set ::ufv::a_function $t
}

# VMD menu registration (when loaded as an extension).
if {![catch {package present vmd}]} {
    catch { vmd_install_extension ufv ufv_gui "Analysis/3D Feature Viewer for UniProt" }
}

puts "\[UFV\] 3D Feature Viewer for UniProt (VMD) loaded - type 'ufv_gui' to open the panel."
puts "      commands: ufv_load, ufv_fetch, ufv_map, ufv_chain, ufv_ptms, ufv_variants, ufv_sites,"
puts "      ufv_mutagenesis, ufv_domains, ufv_topology, ufv_alphamissense, ufv_burden, ufv_plddt,"
puts "      ufv_bfactor, ufv_residue <pos>, ufv_load_structure <key>, ufv_clear, ufv_gui"
