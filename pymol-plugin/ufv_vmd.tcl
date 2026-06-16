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
    # tier cache: recomputed when (uid, molid) pair changes
    variable tiers_uid    ""
    variable tiers_molid  ""
    # persistent rep indices (-1 = no rep currently shown)
    variable lig_rep      -1
    variable focus_rep    -1
    # GUI state
    variable gui_filt_path 1
    variable gui_filt_del  1
    variable gui_filt_ben  0
    variable gui_filt_unc  0
    variable gui_ann_kind  "PTMs"
    variable gui_ann_filter ""
    variable gui_ann_positions {}
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

# ---- reverse numbering (resi → UniProt) -----------------------------------------------------
proc ::ufv::resi_to_uni {chain resi} {
    switch $::ufv::mode {
        identity { return $resi }
        manual {
            variable manual
            if {![info exists manual($chain)]} { return "" }
            lassign $manual($chain) us rs re
            if {$resi < $rs} { return "" }
            if {$re ne "" && $resi > $re} { return "" }
            return [expr {$us + ($resi - $rs)}]
        }
        sifts {
            set var ::ufv::s_$chain
            if {![info exists $var]} { return "" }
            foreach seg [set $var] {
                lassign $seg us ue ps pe
                if {$resi >= $ps && $resi <= $pe} { return [expr {$us + ($resi - $ps)}] }
            }
            return ""
        }
    }
    return $resi
}

# ---- portable temp-file path ----------------------------------------------------------------
proc ::ufv::tempfile {} {
    foreach v {TMPDIR TMP TEMP} {
        if {[info exists ::env($v)] && $::env($v) ne ""} {
            return [file join $::env($v) ufv_[clock clicks].json]
        }
    }
    return /tmp/ufv_[clock clicks].json
}

# ---- export Cα coords as JSON (for --compute-tiers) ----------------------------------------
proc ::ufv::export_ca_coords {molid tmpfile} {
    set sel [atomselect $molid "protein and name CA"]
    set chains [$sel get chain]
    set resids [$sel get resid]
    set xs     [$sel get x]
    set ys     [$sel get y]
    set zs     [$sel get z]
    $sel delete
    set items {}
    foreach ch $chains ri $resids x $xs y $ys z $zs {
        set uni [::ufv::resi_to_uni $ch $ri]
        if {$uni eq ""} continue
        lappend items [format {{"uniPos":%d,"chain":"%s","resi":%d,"x":%.3f,"y":%.3f,"z":%.3f}} \
            $uni $ch $ri $x $y $z]
    }
    set fh [open $tmpfile w]
    puts $fh "\[[join $items ",\n"]\]"
    close $fh
}

# ---- zoom view to a selection ---------------------------------------------------------------
proc ::ufv::zoom_to {molid selstr} {
    set sel [atomselect $molid $selstr]
    if {[$sel num] == 0} { $sel delete; return }
    set center [measure center $sel weight mass]
    $sel delete
    display resetview
    molinfo $molid set center_matrix [transoffset [vecscale -1.0 $center]]
}

# ---- ensure hotspot/hub/pocket tiers are computed and cached --------------------------------
proc ::ufv::ensure_tiers {molid uid} {
    if {$::ufv::tiers_uid eq $uid && $::ufv::tiers_molid eq $molid} { return 1 }
    ::ufv::status "computing hotspots / hubs / pockets..."
    set tmpfile [::ufv::tempfile]
    if {[catch {::ufv::export_ca_coords $molid $tmpfile} err]} {
        ::ufv::status "coord export failed: $err"; return 0
    }
    if {[catch {eval [::ufv::run --compute-tiers $uid $tmpfile]} err]} {
        ::ufv::status "tier computation failed: $err"
        catch { file delete $tmpfile }
        return 0
    }
    catch { file delete $tmpfile }
    set ::ufv::tiers_uid   $uid
    set ::ufv::tiers_molid $molid
    return 1
}

# ---- populate annotation listbox ------------------------------------------------------------
proc ::ufv::refresh_ann_list {} {
    set lb .ufv.canvframe.c.f.ann.lb
    if {![winfo exists $lb]} return
    $lb delete 0 end
    set filter [string tolower $::ufv::gui_ann_filter]
    set positions {}
    switch $::ufv::gui_ann_kind {
        PTMs {
            if {[info exists ::ufv::a_ptms_rich]} {
                foreach e $::ufv::a_ptms_rich {
                    lassign $e pos color desc
                    set line "[format {%5d} $pos]  $desc"
                    if {$filter eq "" || [string first $filter [string tolower $line]] >= 0} {
                        $lb insert end $line
                        lappend positions $pos
                    }
                }
            }
        }
        Variants {
            if {[info exists ::ufv::a_variants_rich]} {
                foreach e $::ufv::a_variants_rich {
                    lassign $e pos color tok wtmut label
                    set skip 0
                    switch $tok {
                        pathogenic  { if {!$::ufv::gui_filt_path} { set skip 1 } }
                        deleterious { if {!$::ufv::gui_filt_del}  { set skip 1 } }
                        benign      { if {!$::ufv::gui_filt_ben}   { set skip 1 } }
                        uncertain   { if {!$::ufv::gui_filt_unc}   { set skip 1 } }
                    }
                    if {$skip} continue
                    set line "[format {%5d} $pos]  $label"
                    if {$filter eq "" || [string first $filter [string tolower $line]] >= 0} {
                        $lb insert end $line
                        lappend positions $pos
                    }
                }
            }
        }
        Sites {
            if {[info exists ::ufv::a_sites_rich]} {
                foreach e $::ufv::a_sites_rich {
                    lassign $e pos color desc
                    set line "[format {%5d} $pos]  $desc"
                    if {$filter eq "" || [string first $filter [string tolower $line]] >= 0} {
                        $lb insert end $line
                        lappend positions $pos
                    }
                }
            }
        }
    }
    set ::ufv::gui_ann_positions $positions
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

# ---- ligands --------------------------------------------------------------------------------
proc ufv_ligands {{molid ""}} {
    set molid [::ufv::resolve $molid]
    if {$::ufv::lig_rep >= 0} {
        catch { mol delrep $::ufv::lig_rep $molid }
        set ::ufv::lig_rep -1
    }
    set chk [atomselect $molid "not protein and not water and noh and not name 'NA' 'CA' 'K' 'MG' 'ZN' 'FE' 'CU' 'MN'"]
    set n [$chk num]; $chk delete
    if {$n == 0} { ::ufv::status "no ligands in mol $molid."; return }
    ::ufv::addrep $molid "Licorice 0.3 12" \
        "not protein and not water and noh and not name NA CA K MG ZN FE CU MN" \
        [::ufv::next_color "#f9a825"]
    set ::ufv::lig_rep [expr {[molinfo $molid get numreps] - 1}]
    ::ufv::status "ligands shown ($n atoms)."
}

proc ufv_ligands_hide {{molid ""}} {
    set molid [::ufv::resolve $molid]
    if {$::ufv::lig_rep >= 0} {
        catch { mol delrep $::ufv::lig_rep $molid }
        set ::ufv::lig_rep -1
    }
    ::ufv::status "ligands hidden."
}

# ---- focus / reset view ---------------------------------------------------------------------
proc ufv_focus {pos {molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    if {$::ufv::focus_rep >= 0} {
        catch { mol delrep $::ufv::focus_rep $molid }
        set ::ufv::focus_rep -1
    }
    set sel [::ufv::sel_positions $molid [list $pos] 0]
    if {$sel eq ""} { ::ufv::status "position $pos not mapped."; return }
    ::ufv::addrep $molid "Licorice 0.3 12" "byres ($sel or within 5 of ($sel))" \
        [::ufv::next_color "#ffa726"]
    set ::ufv::focus_rep [expr {[molinfo $molid get numreps] - 1}]
    ::ufv::zoom_to $molid $sel
    ::ufv::status "focused on UniProt $pos."
}

proc ufv_resetview {{molid ""}} {
    set molid [::ufv::resolve $molid]
    if {$::ufv::focus_rep >= 0} {
        catch { mol delrep $::ufv::focus_rep $molid }
        set ::ufv::focus_rep -1
    }
    set full [atomselect $molid "protein"]
    if {[$full num] > 0} {
        set center [measure center $full weight mass]
        $full delete
        display resetview
        molinfo $molid set center_matrix [transoffset [vecscale -1.0 $center]]
    } else { $full delete; display resetview }
    ::ufv::status "view reset."
}

# ---- align (RMSD fit or centroid fallback) --------------------------------------------------
proc ufv_align {{refmol ""}} {
    set refmol [::ufv::resolve $refmol]
    set moved 0
    foreach m [molinfo list] {
        if {$m == $refmol} continue
        set refsel [atomselect $refmol "protein and backbone"]
        set movsel [atomselect $m     "protein and backbone"]
        if {[$movsel num] > 0 && [$movsel num] == [$refsel num]} {
            if {[catch {
                set mat [measure fit $movsel $refsel]
                set all [atomselect $m all]; $all move $mat; $all delete
            }]} {
                set rc [measure center $refsel weight mass]
                set mc [measure center $movsel weight mass]
                set all [atomselect $m all]; $all moveby [vecsub $rc $mc]; $all delete
            }
        } elseif {[$movsel num] > 0} {
            set rc [measure center $refsel weight mass]
            set mc [measure center $movsel weight mass]
            set all [atomselect $m all]; $all moveby [vecsub $rc $mc]; $all delete
        }
        $refsel delete; $movsel delete
        incr moved
    }
    ::ufv::status "aligned $moved mol(s) to mol $refmol."
}

# ---- geometry-based cartoon colourings ------------------------------------------------------
proc ::ufv::apply_tier_colouring {molid tiervar colors label} {
    if {![info exists $tiervar] || [llength [set $tiervar]] == 0} {
        ::ufv::status "no $label data."; return
    }
    array set buckets {}
    foreach e [set $tiervar] {
        lassign $e pos tier
        lappend buckets($tier) $pos
    }
    foreach {tier color} $colors {
        if {![info exists buckets($tier)]} continue
        set sel [::ufv::sel_positions $molid $buckets($tier) 0]
        if {$sel ne ""} { ::ufv::addrep $molid "NewCartoon" $sel [::ufv::next_color $color] }
    }
    array unset buckets
    ::ufv::status "$label shown ([llength [set $tiervar]] residues)."
}

proc ufv_hotspots {{molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    if {![::ufv::ensure_tiers $molid $::ufv::uid]} return
    ::ufv::apply_tier_colouring $molid ::ufv::a_hotspots \
        {1 #b71c1c 2 #e65100 3 #f9a825} "variant hotspots"
}

proc ufv_contacthubs {{molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    if {![::ufv::ensure_tiers $molid $::ufv::uid]} return
    ::ufv::apply_tier_colouring $molid ::ufv::a_hubs \
        {1 #1a237e 2 #3949ab 3 #7986cb} "contact hubs"
}

proc ufv_pockets {{molid ""}} {
    set molid [::ufv::resolve $molid]; ::ufv::need ""
    if {![::ufv::ensure_tiers $molid $::ufv::uid]} return
    ::ufv::apply_tier_colouring $molid ::ufv::a_pockets \
        {1 #1b5e20 2 #388e3c 3 #81c784} "constraint pockets"
}

proc ufv_show_pocket {positions {molid ""}} {
    set molid [::ufv::resolve $molid]
    set sel [::ufv::sel_positions $molid $positions 0]
    if {$sel eq ""} { ::ufv::status "no pocket residues mapped."; return }
    ::ufv::addrep $molid "Licorice 0.3 12" "byres ($sel)" [::ufv::next_color "#e65100"]
    ::ufv::status "pocket lining shown ([llength $positions] residues)."
}

# ---- CSV export -----------------------------------------------------------------------------
proc ufv_csv {{outfile ""} {molid ""}} {
    ::ufv::need ""
    if {$outfile eq ""} { set outfile "[file normalize ${::ufv::uid}_residue_annotations.csv]" }
    ::ufv::status "writing CSV to $outfile ..."
    if {[catch { eval [::ufv::run --csv $::ufv::uid $outfile] } err]} {
        ::ufv::status "CSV export failed: $err"; return
    }
    ::ufv::status "CSV written: $outfile"
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
    set _reptxt .ufv.canvframe.c.f.rep.t
    if {[winfo exists $_reptxt]} {
        $_reptxt configure -state normal
        $_reptxt delete 1.0 end
        $_reptxt insert end $txt
        $_reptxt configure -state disabled
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

    # Scrollable main area
    frame .ufv.canvframe; pack .ufv.canvframe -fill both -expand 1
    canvas .ufv.canvframe.c -yscrollcommand ".ufv.canvframe.sb set" -highlightthickness 0
    scrollbar .ufv.canvframe.sb -command ".ufv.canvframe.c yview"
    pack .ufv.canvframe.sb -side right -fill y
    pack .ufv.canvframe.c  -side left  -fill both -expand 1
    set f [frame .ufv.canvframe.c.f -padx 10 -pady 8]
    .ufv.canvframe.c create window 0 0 -anchor nw -window $f
    bind $f <Configure> {
        .ufv.canvframe.c configure -scrollregion [.ufv.canvframe.c bbox all]
        .ufv.canvframe.c configure -width [winfo reqwidth .ufv.canvframe.c.f]
    }

    # ---- accession + fetch ---
    set s1 [labelframe $f.s1 -text "Structure" -padx 6 -pady 6]; pack $s1 -fill x
    label $s1.lu -text "UniProt:"; grid $s1.lu -row 0 -column 0 -sticky w
    entry $s1.eu -width 11 -textvariable ::ufv::gui_uid; grid $s1.eu -row 0 -column 1 -sticky w
    set ::ufv::gui_uid $::ufv::uid
    button $s1.fetch -text "Fetch" -command {
        if {$::ufv::gui_uid eq ""} { tk_messageBox -message "Enter an accession."; return }
        ::ufv::status "fetching $::ufv::gui_uid ..."
        ufv_fetch $::ufv::gui_uid
        ::ufv::fetch_structures $::ufv::gui_uid
        .ufv.canvframe.c.f.s1.lb delete 0 end
        foreach st $::ufv::a_structures {
            .ufv.canvframe.c.f.s1.lb insert end "[lindex $st 1]  ([lindex $st 2])"
        }
        ::ufv::ui_function
        ::ufv::refresh_ann_list
        ::ufv::status "fetched $::ufv::gui_uid — choose a structure and Load."
    }
    grid $s1.fetch -row 0 -column 2 -sticky w -padx 4
    button $s1.uniprot -text "Open UniProt ↗" -command {
        if {$::ufv::gui_uid ne ""} {
            catch { exec [auto_execok xdg-open] \
                "https://www.uniprot.org/uniprotkb/$::ufv::gui_uid" & }
            catch { exec [auto_execok open] \
                "https://www.uniprot.org/uniprotkb/$::ufv::gui_uid" & }
        }
    }
    grid $s1.uniprot -row 0 -column 3 -sticky w -padx 2
    listbox $s1.lb -height 5 -width 46 -yscrollcommand "$s1.sb set"
    scrollbar $s1.sb -command "$s1.lb yview"
    grid $s1.lb -row 1 -column 0 -columnspan 4 -sticky ew -pady 4
    grid $s1.sb -row 1 -column 4 -sticky ns
    button $s1.load   -text "Load selected"  -command {
        set i [.ufv.canvframe.c.f.s1.lb curselection]
        if {$i eq ""} { tk_messageBox -message "Select a structure from the list."; return }
        ufv_load_structure [lindex [lindex $::ufv::a_structures $i] 0]
    }
    button $s1.loadall -text "Load all" -command {
        if {![llength $::ufv::a_structures]} { tk_messageBox -message "Fetch an accession first."; return }
        foreach st $::ufv::a_structures { ufv_load_structure [lindex $st 0] }
    }
    button $s1.loadaf -text "Quick AlphaFold" -command {
        if {$::ufv::gui_uid ne ""} { ufv_load $::ufv::gui_uid }
    }
    grid $s1.load    -row 2 -column 0 -columnspan 2 -sticky w -pady 2
    grid $s1.loadall -row 2 -column 2 -sticky w -padx 4 -pady 2
    grid $s1.loadaf  -row 2 -column 3 -sticky w -pady 2

    # ---- numbering ---
    set num [labelframe $f.num -text "Residue numbering" -padx 6 -pady 6]; pack $num -fill x -pady 4
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
            sifts   { ufv_map sifts $::ufv::gui_pdb }
            manual  { ufv_chain $::ufv::gui_ch $::ufv::gui_us $::ufv::gui_rs $::ufv::gui_re }
            default { ufv_map identity }
        }
    }
    grid $num.apply -row 3 -column 0 -columnspan 3 -sticky w -pady 4

    # ---- marker layers ---
    set lay [labelframe $f.lay -text "Layers (markers)" -padx 6 -pady 6]; pack $lay -fill x
    set i 0
    foreach {lbl cmd} {PTMs ufv_ptms "Disease variants" ufv_variants "Func. sites" ufv_sites \
                        Mutagenesis ufv_mutagenesis Ligands ufv_ligands} {
        button $lay.b$i -text $lbl -width 14 \
            -command "if {\$::ufv::gui_uid eq {}} { tk_messageBox -message {Enter an accession.} } \
                      else { ::ufv::need {} ; $cmd }"
        grid $lay.b$i -row [expr {$i/3}] -column [expr {$i%3}] -sticky w -padx 2 -pady 2
        incr i
    }
    # Variant consequence filter
    label $lay.flt -text "Variants filter:"; grid $lay.flt -row 2 -column 0 -sticky w -pady 2
    checkbutton $lay.fp -text "Pathogenic"  -variable ::ufv::gui_filt_path \
        -command { ::ufv::refresh_ann_list }
    checkbutton $lay.fd -text "Deleterious" -variable ::ufv::gui_filt_del \
        -command { ::ufv::refresh_ann_list }
    checkbutton $lay.fb -text "Benign"      -variable ::ufv::gui_filt_ben  \
        -command { ::ufv::refresh_ann_list }
    checkbutton $lay.fu -text "Uncertain"   -variable ::ufv::gui_filt_unc  \
        -command { ::ufv::refresh_ann_list }
    grid $lay.fp -row 2 -column 1 -sticky w
    grid $lay.fd -row 2 -column 2 -sticky w
    grid $lay.fb -row 3 -column 1 -sticky w
    grid $lay.fu -row 3 -column 2 -sticky w

    # ---- cartoon colouring ---
    set col [labelframe $f.col -text "Cartoon colouring" -padx 6 -pady 6]; pack $col -fill x -pady 4
    set i 0
    foreach {lbl cmd} {Domains ufv_domains Topology ufv_topology AlphaMissense ufv_alphamissense \
                        Burden ufv_burden "pLDDT (pred.)" ufv_plddt "B-factor (exp.)" ufv_bfactor \
                        Hotspots ufv_hotspots "Contact hubs" ufv_contacthubs "Constraint pocket" ufv_pockets} {
        button $col.b$i -text $lbl -width 16 \
            -command "if {\$::ufv::gui_uid eq {}} { tk_messageBox -message {Enter an accession.} } \
                      else { ::ufv::need {} ; $cmd }"
        grid $col.b$i -row [expr {$i/3}] -column [expr {$i%3}] -sticky w -padx 2 -pady 2
        incr i
    }

    # ---- residue report ---
    set rep [labelframe $f.rep -text "Residue report" -padx 6 -pady 6]; pack $rep -fill x -pady 4
    label $rep.l -text "UniProt pos:"; grid $rep.l -row 0 -column 0 -sticky w
    entry $rep.e -width 7 -textvariable ::ufv::gui_pos; grid $rep.e -row 0 -column 1 -sticky w
    button $rep.show -text "Show" -command {
        if {$::ufv::gui_pos eq ""} { tk_messageBox -message "Enter a residue position."; return }
        ::ufv::status "fetching residue $::ufv::gui_pos ..."
        ufv_residue $::ufv::gui_pos
        ::ufv::status "residue $::ufv::gui_pos report shown."
    }
    button $rep.focus -text "Focus 3D" -command {
        if {$::ufv::gui_pos eq ""} { tk_messageBox -message "Enter a residue position."; return }
        ufv_focus $::ufv::gui_pos
    }
    button $rep.reset -text "Reset view" -command { ufv_resetview }
    grid $rep.show  -row 0 -column 2 -sticky w -padx 2
    grid $rep.focus -row 0 -column 3 -sticky w -padx 2
    grid $rep.reset -row 0 -column 4 -sticky w -padx 2
    text $rep.t -height 9 -width 52 -wrap word -state disabled -font {Courier 9} \
        -yscrollcommand "$rep.sb set"
    scrollbar $rep.sb -command "$rep.t yview"
    grid $rep.t  -row 1 -column 0 -columnspan 5 -sticky nsew -pady 4
    grid $rep.sb -row 1 -column 5 -sticky ns
    grid rowconfigure $rep 1 -weight 1

    # ---- annotation list ---
    set ann [labelframe $f.ann -text "Annotations" -padx 6 -pady 6]; pack $ann -fill x -pady 4
    frame $ann.top; pack $ann.top -fill x
    foreach {lbl kind} {PTMs PTMs Variants Variants Sites Sites} {
        radiobutton $ann.top.r$kind -text $lbl -variable ::ufv::gui_ann_kind -value $kind \
            -command ::ufv::refresh_ann_list
        pack $ann.top.r$kind -side left
    }
    label $ann.top.fl -text "  Filter:"; pack $ann.top.fl -side left
    entry $ann.top.fe -width 14 -textvariable ::ufv::gui_ann_filter
    pack $ann.top.fe -side left
    bind $ann.top.fe <KeyRelease> { ::ufv::refresh_ann_list }
    listbox $ann.lb -height 8 -width 52 -yscrollcommand "$ann.sb set" -font {Courier 9}
    scrollbar $ann.sb -command "$ann.lb yview"
    pack $ann.lb -side left  -fill both -expand 1 -pady 4
    pack $ann.sb -side right -fill y    -pady 4
    bind $ann.lb <<ListboxSelect>> {
        set idx [.ufv.canvframe.c.f.ann.lb curselection]
        if {$idx ne "" && $idx < [llength $::ufv::gui_ann_positions]} {
            set pos [lindex $::ufv::gui_ann_positions $idx]
            set ::ufv::gui_pos $pos
            ufv_focus $pos
        }
    }

    # ---- function context ---
    label $f.func -textvariable ::ufv::a_function -anchor w -justify left -wraplength 380 -fg "#5b3a8e"
    pack $f.func -fill x

    # ---- bottom buttons + status ---
    label .ufv.status -text "Ready." -anchor w -relief sunken -padx 6
    pack .ufv.status -fill x -side bottom
    frame $f.bot; pack $f.bot -fill x -pady 4
    button $f.bot.align -text "Align mols" -command { ufv_align }
    button $f.bot.clear -text "Clear overlays" -command { ufv_clear }
    button $f.bot.csv   -text "Export CSV" -command {
        if {$::ufv::uid eq ""} { tk_messageBox -message "Fetch an accession first."; return }
        set fn [tk_getSaveFile -defaultextension .csv \
            -initialfile "${::ufv::uid}_residue_annotations.csv" \
            -filetypes {{"CSV files" .csv} {"All files" *}}]
        if {$fn ne ""} { ufv_csv $fn }
    }
    pack $f.bot.align -side left -padx 2
    pack $f.bot.clear -side left -padx 2
    pack $f.bot.csv   -side left -padx 2
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
puts "      ufv_mutagenesis, ufv_ligands, ufv_ligands_hide, ufv_domains, ufv_topology,"
puts "      ufv_alphamissense, ufv_burden, ufv_plddt, ufv_bfactor, ufv_hotspots, ufv_contacthubs,"
puts "      ufv_pockets, ufv_focus <pos>, ufv_resetview, ufv_align, ufv_show_pocket <positions>,"
puts "      ufv_residue <pos>, ufv_csv [outfile], ufv_load_structure <key>, ufv_clear, ufv_gui"
