#!/usr/bin/env python3
import math

def get_mods(kmer_list, d):
    assert int(d) > 0
    mod_list = []
    for kmer in kmer_list:
        if kmer % d == 0:
            mod_list.append(kmer)

    return mod_list

def partition_windows(kmer_list):
    table_size = math.floor(len(kmer_list)/1000)
    thousand_dict = {}
    for i in range(1000):
        start_size = i*table_size
        end_size = (i*table_size + table_size - 1)
        name = f"{start_size}-{end_size}"
        thousand_dict[name] = kmer_list[start_size:end_size]
    return thousand_dict

def binomial_distance(containment_value,kmer_value):
    return math.pow(containment_value, (1/kmer_value))

def containment(set1,set2):
    intersection = set1.intersection(set2)
    try:
        return float(len(intersection)/len(set1))
    except:
        return 0