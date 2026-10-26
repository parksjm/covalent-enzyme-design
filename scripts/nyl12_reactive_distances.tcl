# Native VMD bond labels update their distance values with the current frame.
namespace eval ::nyl12_reactive {
    variable frame_stages
    variable stage_pairs
    variable label_ids
    variable dynamic_reps
    variable current_stage
}
proc ::nyl12_reactive::install {molid stages pairs_by_stage {dynamic_by_stage {}} {persistent_pairs {}} {proton_atom {}} {proton_name {}} {persistent_angles {}}} {
    variable frame_stages
    variable stage_pairs
    variable label_ids
    variable dynamic_reps
    variable current_stage
    if {[info exists frame_stages($molid)]} {error "Reactive labels already installed"}
    if {[llength $stages] != [molinfo $molid get numframes]} {error "Reaction stage/frame mismatch"}
    set frame_stages($molid) $stages
    set stage_pairs($molid) $pairs_by_stage
    set label_ids($molid) [dict create]
    set dynamic_reps($molid) [dict create]
    set current_stage($molid) ""
    dict for {stage pairs} $pairs_by_stage {
        foreach pair $pairs {
            lassign $pair a name_a b name_b
            foreach index [list $a $b] name [list $name_a $name_b] {
                set selection [atomselect $molid "index $index and name $name"]
                set count [$selection num]
                $selection delete
                if {$count != 1} {error "Reactive atom identity mismatch: $index/$name"}
            }
            set key [join [lsort -integer [list $a $b]] ,]
            if {![dict exists $label_ids($molid) $key]} {
                set id [llength [label list Bonds]]
                label add Bonds $molid/$a $molid/$b
                dict set label_ids($molid) $key $id
                label hide Bonds $id
            }
        }
    }
    dict for {stage pairs} $dynamic_by_stage {
        foreach pair $pairs {
            lassign $pair a name_a b name_b cutoff
            foreach index [list $a $b] name [list $name_a $name_b] {
                set selection [atomselect $molid "index $index and name $name"]
                set count [$selection num]
                $selection delete
                if {$count != 1} {error "Dynamic bond atom mismatch: $index/$name"}
            }
            if {$cutoff <= 0} {error "Dynamic-bond cutoff must be positive"}
            set first [atomselect $molid "index $a"]
            set first_bonds [lindex [$first getbonds] 0]
            set first_position [lsearch -exact $first_bonds $b]
            if {$first_position >= 0} {
                $first setbonds [list [lreplace $first_bonds $first_position $first_position]]
            }
            $first delete
            set second [atomselect $molid "index $b"]
            set second_bonds [lindex [$second getbonds] 0]
            set second_position [lsearch -exact $second_bonds $a]
            if {$second_position >= 0} {
                $second setbonds [list [lreplace $second_bonds $second_position $second_position]]
            }
            $second delete
            set rep [mol numreps $molid]
            mol representation DynamicBonds $cutoff 0.12 12
            mol color ColorID 1
            mol selection "index $a $b"
            mol material Opaque
            mol addrep $molid
            mol modselect $rep $molid "index $a $b"
            mol selupdate $rep $molid on
            mol showrep $molid $rep off
            dict lappend dynamic_reps($molid) $stage $rep
        }
    }
    foreach pair $persistent_pairs {
        lassign $pair title a name_a b name_b
        foreach index [list $a $b] name [list $name_a $name_b] {
            set selection [atomselect $molid "index $index and name $name"]
            set count [$selection num]
            $selection delete
            if {$count != 1} {error "Persistent label atom mismatch: $index/$name"}
        }
        set id [llength [label list Bonds]]
        label add Bonds $molid/$a $molid/$b
        label show Bonds $id
    }
    foreach angle $persistent_angles {
        lassign $angle title a name_a b name_b c name_c
        foreach index [list $a $b $c] name [list $name_a $name_b $name_c] {
            set selection [atomselect $molid "index $index and name $name"]
            set count [$selection num]
            $selection delete
            if {$count != 1} {error "Persistent angle atom mismatch: $index/$name"}
        }
        if {[llength [lsort -unique -integer [list $a $b $c]]] != 3} {
            error "Persistent angle needs three distinct atoms"
        }
        set id [llength [label list Angles]]
        label add Angles $molid/$a $molid/$b $molid/$c
        label show Angles $id
    }
    if {$proton_atom ne ""} {
        set proton [atomselect $molid "index $proton_atom and name $proton_name"]
        set count [$proton num]
        $proton delete
        if {$count != 1} {error "Proton sphere atom mismatch: $proton_atom/$proton_name"}
        set rep [mol numreps $molid]
        mol representation VDW 0.32 16
        mol color ColorID 7
        mol selection "index $proton_atom"
        mol material Opaque
        mol addrep $molid
        mol modselect $rep $molid "index $proton_atom"
    }
    trace add variable ::vmd_frame($molid) write [list ::nyl12_reactive::update $molid]
    ::nyl12_reactive::update $molid
}
proc ::nyl12_reactive::update {molid args} {
    variable frame_stages
    variable stage_pairs
    variable label_ids
    variable dynamic_reps
    variable current_stage
    set frame [molinfo $molid get frame]
    set stage [lindex $frame_stages($molid) $frame]
    if {$stage eq ""} {error "No reaction stage for frame $frame"}
    if {$stage eq $current_stage($molid)} {return}
    foreach id [dict values $label_ids($molid)] {label hide Bonds $id}
    foreach pair [dict get $stage_pairs($molid) $stage] {
        lassign $pair a name_a b name_b
        set key [join [lsort -integer [list $a $b]] ,]
        label show Bonds [dict get $label_ids($molid) $key]
    }
    foreach reps [dict values $dynamic_reps($molid)] {
        foreach rep $reps {mol showrep $molid $rep off}
    }
    if {[dict exists $dynamic_reps($molid) $stage]} {
        foreach rep [dict get $dynamic_reps($molid) $stage] {
            mol showrep $molid $rep on
        }
    }
    set current_stage($molid) $stage
}
