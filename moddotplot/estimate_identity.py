#!/usr/bin/env python3
import math

def divide_into_chunks(lst, res):
    chunk_size = math.floor(len(lst) / res)
    chunks = [lst[i:i+chunk_size] for i in range(0, len(lst), chunk_size)]
    while len(chunks) > res:
        chunks.pop()
    while len(chunks) < res:
        chunks.append([])
    return chunks

def get_mods(kmer_list, s, res):
    assert int(s) > 0
    mod_list = []
    mod_list_prep = divide_into_chunks(kmer_list, res)
    for element in mod_list_prep:
        mod_list.append([])
        for kmer in element:
            if kmer % s == 0:
                mod_list[-1].append(kmer)
    return mod_list

def partition_evenly_spaced_modimizers(mod_list, seq_length, resolution):
    thousand_dict = {}
    for i in range(resolution):
        start_site = i * round(seq_length/resolution)
        end_site = (i * round(seq_length/resolution) + round(seq_length/resolution) - 1)
        name = f"{start_site}-{end_site}"
        thousand_dict[name] = mod_list[i]
    return thousand_dict

def partition_windows(kmer_list, seq_length, resolution, sparsity):
    table_size = round(len(kmer_list)/resolution)
    actual_ratio = seq_length/len(kmer_list)
    thousand_dict = {}
    for i in range(resolution):
        start_size = i*table_size
        end_size = (i*table_size + table_size - 1)
        start_site = i * round(seq_length/resolution)
        end_site = (i * round(seq_length/resolution) + round(seq_length/resolution) - 1)
        name = f"{start_site}-{end_site}"
        thousand_dict[name] = kmer_list[start_size:end_size]
    return thousand_dict


def partition_windows_limit(kmer_list, seq_length, resolution, upper_limit, real_limit):
    table_size = round(len(kmer_list)/resolution)
    print("-----")
    threshold_amount = round(real_limit/upper_limit * 1000)
    print(threshold_amount)
    print("-----")
    thousand_dict = {}
    for i in range(threshold_amount):
        start_size = i*table_size
        end_size = (i*table_size + table_size - 1)
        start_site = i * round(seq_length/threshold_amount)
        end_site = (i * round(seq_length/threshold_amount) + round(seq_length/threshold_amount) - 1)
        name = f"{start_site}-{end_site}"
        thousand_dict[name] = kmer_list[start_size:end_size]
    if threshold_amount < resolution:
        for i in range(threshold_amount, resolution):
            start_size = i*table_size
            end_size = (i*table_size + table_size - 1)
            start_site = i * round(seq_length/threshold_amount)
            end_site = (i * round(seq_length/threshold_amount) + round(seq_length/threshold_amount) - 1)
            name = f"{start_site}-{end_site}"
            thousand_dict[name] = []

    return thousand_dict

def create_coordinates(kmer_list, resolution):
    zoomed_limit = math.floor(len(kmer_list)/resolution)
    tmp = []
    for k in range(resolution):
        tmp.append(set())
        for j in range(zoomed_limit):
            tmp[k].add(kmer_list[(k*zoomed_limit)+j])

    return tmp

def create_coordinates2(kmer_list):
    tmp = []
    for k in kmer_list:
        tmp.append(set())
        for j in k:
            tmp[-1].add(j)

    return tmp

def binomial_distance(containment_value,kmer_value):
    return math.pow(containment_value, (1/kmer_value))

def containment(set1,set2):
    intersection = set1.intersection(set2)
    try:
        if len(set1) > len(set2):
            return float(len(intersection)/len(set1))
        else:
            return float(len(intersection)/len(set2))
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