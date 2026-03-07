#!/bin/bash

cat $1 | awk -F '\t' '$3=="gene" {
    id = name = ""
    if (match($9, /gene_id "[^"]+"/)) {
        id = substr($9, RSTART+9, RLENGTH-10)
    }
    if (match($9, /gene_name "[^"]+"/)) {
        name = substr($9, RSTART+11, RLENGTH-12)
    }
    if (id != "" && name != "") {
        print id "\t" name
    }
}' | sort -u
