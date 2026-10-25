# Native VMD bond labels update their distance values with the current frame.
namespace eval ::nyl12_reactive {
    variable frame_stages
    variable stage_pairs
    variable label_ids
    variable current_stage
}
proc ::nyl12_reactive::install {molid stages pairs_by_stage} {
    variable frame_stages
    variable stage_pairs
    variable label_ids
    variable current_stage
    if {[info exists frame_stages($molid)]} {error "Reactive labels already installed"}
    if {[llength $stages] != [molinfo $molid get numframes]} {error "Reaction stage/frame mismatch"}
    set frame_stages($molid) $stages
    set stage_pairs($molid) $pairs_by_stage
    set label_ids($molid) [dict create]
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
    trace add variable ::vmd_frame($molid) write [list ::nyl12_reactive::update $molid]
    ::nyl12_reactive::update $molid
}
proc ::nyl12_reactive::update {molid args} {
    variable frame_stages
    variable stage_pairs
    variable label_ids
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
    set current_stage($molid) $stage
}
