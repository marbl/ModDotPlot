#!/usr/bin/env python3
import math

def get_mods(kmer_list, d):
    assert int(d) > 0
    mod_list = []
    for kmer in kmer_list:
        if kmer % d == 0:
            mod_list.append(kmer)

    return mod_list

def partition_windows(kmer_list, seq_length, resolution):
    table_size = math.floor(len(kmer_list)/resolution)
    thousand_dict = {}
    for i in range(resolution):
        start_size = i*table_size
        end_size = (i*table_size + table_size - 1)
        start_site = i * round(seq_length/resolution)
        end_site = (i * round(seq_length/resolution) + round(seq_length/resolution) - 1)
        name = f"{start_site}-{end_site}"
        thousand_dict[name] = kmer_list[start_size:end_size]
    return thousand_dict

def create_coordinates(kmer_list, resolution):
    zoomed_limit = math.floor(len(kmer_list)/resolution)
    tmp = []
    for k in range(resolution):
        tmp.append(set())
        for j in range(zoomed_limit):
            tmp[k].add(kmer_list[(k*zoomed_limit)+j])

    return tmp

def binomial_distance(containment_value,kmer_value):
    return math.pow(containment_value, (1/kmer_value))

def containment(set1,set2):
    intersection = set1.intersection(set2)
    try:
        return float(len(intersection)/len(set1))
    except:
        return 0

def jaccard(set1,set2):
    intersection = set1.intersection(set2)
    try:
        return float(len(intersection)/(len(set1) + len(set2) - len(intersection)))
    except:
        return 0

def poisson_distance(jaccard_value,kmer_value):
    return (-1/kmer_value)*math.log((2*jaccard_value/(1+jaccard_value)))