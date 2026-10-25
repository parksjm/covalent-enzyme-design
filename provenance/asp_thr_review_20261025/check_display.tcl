# Mock VMD command dispatch only; real DCD/topology data are checked separately.
set package [file normalize [lindex $argv 0]]
set work [file normalize [lindex $argv 1]]
file mkdir $work
foreach name {nyl12_es_optimization.vmd nyl12_es_initial_final.vmd nyl12_es_vmd_common.tcl qm_indices.tcl} {
    if {$name eq "nyl12_es_vmd_common.tcl"} {file copy [lindex $argv 2] $work/$name} else {file copy $package/$name $work/$name}
}
foreach name {stripped.parm7.gz md_solvated.parm7.gz optimization.dcd initial_final.dcd} {
    set f [open $work/$name w]; puts $f "mock input"; close $f
}
rename exec actual_exec
proc exec {args} {set f [open [lindex $args end] w]; puts $f "mock topology"; close $f}
proc mol {action args} {
    lappend ::mol_calls [linsert $args 0 $action]
    if {$action eq "new"} {
        if {[file tail [lindex $args 0]] eq "stripped.parm7"} {set ::mock_atoms $::nyl12_stripped_atom_count} else {set ::mock_atoms 103251}
        return 0
    }
    if {$action eq "addfile"} {
        if {[file tail [lindex $args 0]] eq "optimization.dcd"} {set ::mock_frames 48} else {set ::mock_frames 2}
    }
}
proc molinfo {id operation args} {
    if {$operation eq "get"} {
        if {[lindex $args 0] eq "numatoms"} {return $::mock_atoms}
        if {[lindex $args 0] eq "numframes"} {return $::mock_frames}
    }
}
proc atomselect {args} {return mock_selection}
proc mock_selection {action args} {if {$action eq "num"} {return 98}}
proc measure {args} {return {0 0 0}}
proc trans {args} {return {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}
proc display {args} {}
proc scale {args} {}
proc animate {args} {}
source $work/nyl12_es_optimization.vmd
if {$::mock_frames != 48 || $::mock_atoms != $::nyl12_stripped_atom_count} {error "default movie mismatch"}
initial; final; pocket; system
full
if {$::mock_frames != 2 || $::mock_atoms != 103251} {error "full comparison mismatch"}
initial; final
file mkdir $work/missing
file copy $package/nyl12_es_optimization.vmd $work/missing/nyl12_es_optimization.vmd
if {![catch {source $work/missing/nyl12_es_optimization.vmd} message] || ![string match "Missing viewer input:*" $message]} {error "Missing input did not fail explicitly"}
puts "Mocked Tcl loader/helper and missing-input checks passed; no native VMD rendering tested."

if {[lsearch -exact $::mol_calls {showrep 0 0 off}] < 0} {error "Protein cartoon not hidden"}
if {[lsearch -exact $::mol_calls {color ColorID 3}] < 0} {error "Orange overlay missing"}
set found 0
foreach call $::mol_calls {if {[lindex $call 0] eq "selection" && [string match {resname LIG and element C and not (*)} [lindex $call 1]]} {set found 1}}
if {!$found} {error "Non-QM carbon selection missing"}
puts "Orange non-QM carbon overlay and hidden cartoon checks passed."
