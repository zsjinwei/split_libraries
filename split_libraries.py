# -*- coding: utf-8 -*-
# author： huangjinwei & liaoying
# 2019-03-10

import subprocess
import os, sys, re
import linecache
import datetime
from tqdm import tqdm

# 删除非空文件夹
def delete_dire(dire):
    dir_list = []
    for root, dirs, files in os.walk(dire):
        for afile in files:
            os.remove(os.path.join(root, afile))
        for adir in dirs:
            dir_list.append(os.path.join(root, adir))
    for bdir in dir_list:
        os.rmdir(bdir)

# input_file: 输入barcode fastq文件
# barcodes_need_count: barcode类别数量
# barcode_len: barcode类别长度
def extract_barcode(input_file, barcodes_need_count, barcode_len, barcode_line):
    count = 0
    barcodes = {}
    item = []
    fd = open(input_file, 'r')
    for line in tqdm(fd.readlines(), desc="extract_barcode " + input_file):
        # print(line.strip())
        item.append(line)
        if (count+1) % barcode_line == 0:
            barcode = item[1][0:barcode_len]

            p = re.compile(r'^@.*?:.*?:.*?:.*?:(.*?:.*?:.*?) ')
            barcode_id = p.findall(item[0])[0]
            if(barcode_id == ""):
                print("Got empty barcode_id")
                item = []
                continue
            if barcode in barcodes.keys():
                barcodes[barcode].append(barcode_id + "*" + str(count - 3))
            else:
                barcodes[barcode] = []
                barcodes[barcode].append(barcode_id + "*" + str(count - 3))
            item = []

        count += 1

    barcodes_sorted = sorted(barcodes.items(), key=lambda d:len(d[1]), reverse = True)
    barcodes_need = barcodes_sorted[0:barcodes_need_count]

    i = 0
    correct_count = 0
    error_count = 0
    error_barcode_id = []
    for it in barcodes_sorted:
        if(i < barcodes_need_count):
            correct_count += len(it[1])
        else:
            error_count += len(it[1])
            for id_line_it in it[1]:
                id_line = id_line_it.split('*')
                id = id_line[0]
                line = id_line[1]
                error_barcode_id.append(id)
        i += 1

    # print(error_barcode_id)

    fd.close()

    print("Correct rate = " + str((correct_count / (correct_count + error_count))*100) + "% (" + str(correct_count) + "/" + str(error_count) + ")")

    # print(barcodes_need)
    return barcodes_need, error_barcode_id

# type:id表转id_type表
def id_match_type(type_match_id):
    #id_fd = open('id_dup.txt', 'w+')
    id_match_type = {}
    idcounter = 0
    for item in tqdm(type_match_id, desc="id_match_type"):
        _type = item[0]
        for i in range(len(item[1])):
            id_line = item[1][i].split('*')
            id = id_line[0]
            line = id_line[1]
            if id in id_match_type.keys():
                #re = linecache.getlines(input_file)[id_match_type[id][1]:id_match_type[id][1] + 4]
                #for it in re:
                #    id_fd.write(it)
                #re = linecache.getlines(input_file)[int(line) : int(line) + 4]
                #for it in re:
                #    id_fd.write(it)
                idcounter += 1
            id_match_type[id] = [_type, int(line)]
    print("duplicate id count = " + str(idcounter))
    #id_fd.close()
    return id_match_type

# barcodes_R1: 从文件1(R1或R2)得到的type-id表
# id_match_type_R2: 从文件2得到的type-id转换后的id-type表
# input_file: 文件1输入文件（将会从此文件读取指定行输出到结果文件）
# output_path: 输出文件路径
def split_to_type_file(barcodes_R1, id_match_type_R2, barcode_id_union, file_sample_meta_data, input_file, output_path, out_suffix, barcode_line):
    found_id_count = 0
    miss_id_count = 0
    open_files = {}

    id_match_type_R2_set = set(id_match_type_R2.keys())
    for item in tqdm(barcodes_R1, desc="split_to_type_file " + input_file):
        R1_type = item[0]
        R1_ids = item[1]
        R1_ids.sort()
        for i in range(len(R1_ids)):
            R1_id_line = item[1][i].split('*')
            R1_id = R1_id_line[0]
            R1_line = int(R1_id_line[1])

            #if R1_id not in barcode_id_union:
            #    continue
            #else:
                # print(R1_id + " is error barcode")
            #    pass

            #sprint("find barcode_id_union spend " + str((endtime - starttime).seconds) + " sec")
            if R1_id in id_match_type_R2_set:
                R2_type_line = id_match_type_R2[R1_id]
                R2_type = R2_type_line[0]
                R2_line = R2_type_line[1]
                re = linecache.getlines(input_file)[R1_line : R1_line + barcode_line]
                # file_path = output_path + '/' + R1_type + "_" + R2_type + ".fastq"
                if (R1_type + R2_type) not in file_sample_meta_data.keys():
                    print((R1_type + R2_type) + " is not in metadata")
                    continue
                file_path = output_path + '/' + file_sample_meta_data[R1_type + R2_type] + out_suffix + ".fastq"
                if file_path not in open_files:
                    wr_fd = open(file_path, "w+")
                    open_files[file_path] = wr_fd
                else:
                    wr_fd = open_files[file_path]

                #wr_fd.write(id + " " + R1_type + " " + str(R2_line) + "\n")
                if R1_id in barcode_id_union:
                    for it in re:
                        wr_fd.write(it)

                found_id_count += 1
            else:
                # print("Warning: id = " + R1_id + " is not found")
                miss_id_count += 1
                pass

    for key in open_files:
        open_files[key].close()

    return found_id_count, miss_id_count

def load_metadata(meta_file, pass_line, barcode_len):
    line_count = 0
    type_sample = {}
    fr = open(meta_file, 'r')
    for line in tqdm(fr.readlines(), desc="load metadata"):
        line_count += 1
        if line_count <= pass_line:
            continue
        lineArr = line.strip().split()
        type_sample[lineArr[1]] = lineArr[0]
        type_sample[lineArr[1][barcode_len:] + lineArr[1][:barcode_len]] = lineArr[0]
    fr.close()
    return type_sample

def get_error_barcode(id_match_type_R1, id_match_type_R2):
    id_match_type_R1_set = set(id_match_type_R1.keys())
    id_match_type_R2_set = set(id_match_type_R2.keys())
    return set(id_match_type_R1_set).intersection(set(id_match_type_R2_set))


def usage(script_name):
    # python split_libraries.py L1P1_R1.fastq L1P1_R2.fastq metadata-L1P1.txt ./result 4 12 12 8 2 _R1 _R2
    print("python " + script_name +" R1_file_path R2_file_path Meta_file_path output_path barcode_line barcode_len R1_type_count R2_type_count meta_pass_line R1_output_suffix R2_output_suffix")

if __name__ == "__main__":

    # R1_file = 'L1P1_R1.fastq' # 输入R1文件名
    # R2_file = 'L1P1_R2.fastq' # 输入R2文件名
    # metadata_file = 'metadata-L1P1.txt' # 输入metadata文件名
    # output_path = './result' # 输出路径
    # barcode_line = 4 # 每个barcode占多少行
    # barcode_len = 12 # barcode长度
    # R1_type_count = 12 # R1文件包含barcode类型数量
    # R2_type_count = 8 # R2文件包含barcode类型数量
    # meta_pass_line = 2 # meta文件里面多少行文件头
    # R1_output_suffix = "_R1" # 输出R1文件名后缀
    # R2_output_suffix = "_R2" # 输出R2文件名后缀

    if(len(sys.argv) != 12):
        usage(sys.argv[0])
        exit()

    R1_file = sys.argv[1]
    R2_file = sys.argv[2]
    metadata_file = sys.argv[3]
    output_path = sys.argv[4]
    barcode_line = int(sys.argv[5])
    barcode_len = int(sys.argv[6])
    R1_type_count = int(sys.argv[7])
    R2_type_count = int(sys.argv[8])
    meta_pass_line = int(sys.argv[9])
    R1_output_suffix = sys.argv[10]
    R2_output_suffix = sys.argv[11]

    print("R1_file = " + str(R1_file))
    print("R2_file = " + str(R2_file))
    print("metadata_file = " + str(metadata_file))
    print("output_path = " + str(output_path))
    print("barcode_line = " + str(barcode_line))
    print("barcode_len = " + str(barcode_len))
    print("R1_type_count = " + str(R1_type_count))
    print("R2_type_count = " + str(R2_type_count))
    print("meta_pass_line = " + str(meta_pass_line))
    print("R1_output_suffix  = " + str(R1_output_suffix))
    print("R2_output_suffix  = " + str(R2_output_suffix))

    # R1_file = 'part_R1_part.fastq'
    # R2_file = 'part_R2_part.fastq'

    starttime = datetime.datetime.now()
    global_starttime = starttime

    # 提取barcode类型
    print("\n\nbegin to extract barcode")
    barcodes_R1_need, R1_error_barcode_id = extract_barcode(R1_file, R1_type_count, barcode_len, barcode_line)
    barcodes_R2_need, R2_error_barcode_id = extract_barcode(R2_file, R2_type_count, barcode_len, barcode_line)

    endtime = datetime.datetime.now()
    print("extract barcode spend " + str((endtime - starttime).seconds) + " sec")

    starttime = datetime.datetime.now()

    # 求并集
    # print("\n\nbegin to merge error barcode id")
    # barcode_id_union = set(R1_error_barcode_id).union(set(R2_error_barcode_id))
    # print(barcode_id_union)

    endtime = datetime.datetime.now()
    print("merge error barcode spend " + str((endtime - starttime).seconds) + " sec")

    starttime = datetime.datetime.now()

    print("\n\nbegin to trans type-id table")
    id_match_type_R1 = id_match_type(barcodes_R1_need)
    id_match_type_R2 = id_match_type(barcodes_R2_need)

    endtime = datetime.datetime.now()
    print("trans type-id table spend " + str((endtime - starttime).seconds) + " sec")

    barcode_id_union = get_error_barcode(id_match_type_R1, id_match_type_R2)

    starttime = datetime.datetime.now()

    print("\n\nbegin to load metadata")
    type_sample = load_metadata(metadata_file, meta_pass_line, barcode_len)

    endtime = datetime.datetime.now()
    print("slpit to type file spend " + str((endtime - starttime).seconds) + " sec")

    starttime = datetime.datetime.now()

    print("\n\nbegin to split to type file")
    delete_dire(output_path)

    os.makedirs(output_path + "/R1_R2")
    os.makedirs(output_path + "/R2_R1")

    R1_found_id_count, R1_miss_id_count = split_to_type_file(barcodes_R1_need, id_match_type_R2, barcode_id_union, type_sample, R1_file, output_path + "/R1_R2/", R1_output_suffix, barcode_line)
    R2_found_id_count, R2_miss_id_count = split_to_type_file(barcodes_R2_need, id_match_type_R1, barcode_id_union, type_sample, R2_file, output_path + "/R2_R1/", R2_output_suffix, barcode_line)

    endtime = datetime.datetime.now()

    print("slpit to type file spend " + str((endtime - starttime).seconds) + " sec")
    print("\n\nTotal time = " + str((endtime - global_starttime).seconds) + " sec")
    print("\nGot R1 = " + str(R1_found_id_count) + " found | " + str(R1_miss_id_count) + " miss and R2 = " + str(R2_found_id_count) + " found | " + str(R2_miss_id_count) + " miss")
    print("\nDone!")
    exit()