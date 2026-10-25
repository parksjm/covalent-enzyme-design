proc nyl12_load {topology filename expected_frames indices expected_atoms} {
    set top [file join $::nyl12_dir $topology]
    foreach path [list "$top.gz" [file join $::nyl12_dir $filename]] {
        if {![file isfile $path]} {error "Missing VMD input: $path"}
    }
    if {![file isfile $top]} {
        if {[auto_execok gzip] eq ""} {error "gzip is required to expand the topology"}
        exec gzip -dc "$top.gz" > $top
    }
    set ::nyl12_mol [mol new $top type parm7 waitfor all]
    mol addfile [file join $::nyl12_dir $filename] type dcd waitfor all molid $::nyl12_mol
    if {[molinfo $::nyl12_mol get numatoms] != $expected_atoms || [molinfo $::nyl12_mol get numframes] != $expected_frames} {
        error "Unexpected full-system atom/frame count"
    }
    set ::nyl12_current_indices $indices
    set qm "index [join $indices { }]"
    set check [atomselect $::nyl12_mol $qm]
    if {[$check num] != 98} {error "Expected all 98 physical QM atoms"}
    $check delete
    mol delrep 0 $::nyl12_mol
    foreach selection [list protein {resname LIG} $qm "same residue as (water and within 5 of ($qm))" {water}] style {NewCartoon Licorice Licorice Licorice Lines} {
        mol representation $style
        mol color Element
        mol selection $selection
        mol material Opaque
        mol addrep $::nyl12_mol
    }
    # Keep complete ligand bonds underneath a slightly wider carbon overlay.
    mol representation Licorice 0.31 12.0 12.0
    mol color ColorID 3
    mol selection "resname LIG and element C and not ($qm)"
    mol material Opaque
    mol addrep $::nyl12_mol
    mol showrep $::nyl12_mol 0 off
    mol selupdate 3 $::nyl12_mol on
    mol showrep $::nyl12_mol 4 off
    mol rename $::nyl12_mol "Nyl12 PA6: $filename"
    display depthcue off
    display resetview
    pocket
}
proc initial {} {animate goto 0 $::nyl12_mol}
proc final {} {animate goto [expr {[molinfo $::nyl12_mol get numframes] - 1}] $::nyl12_mol}
proc pocket {} {
    mol showrep $::nyl12_mol 4 off
    set view [atomselect $::nyl12_mol "resname LIG or index [join $::nyl12_current_indices { }]"]
    display resetview
    molinfo $::nyl12_mol set center_matrix [trans origin [measure center $view]]
    $view delete
    scale by 3.0
}
proc system {} {mol showrep $::nyl12_mol 4 on; display resetview}
proc optimization {} {nyl12_load stripped.parm7 optimization.dcd 48 $::nyl12_qm_indices $::nyl12_stripped_atom_count}
proc full {} {nyl12_load md_solvated.parm7 initial_final.dcd 2 $::nyl12_full_qm_indices 103251}
