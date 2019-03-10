# -*- coding: utf-8 -*-
# author： huangjinwei & liaoying
# 2019-03-10

import subprocess
import os, sys, re
import linecache
import datetime

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
    for line in fd.readlines():
        # print(line.strip())
        item.append(line)
        if (count+1) % barcode_line == 0:
            barcode = item[1][0:barcode_len]

            p = re.compile(r'^@.*?:.*?:.*?:.*?:.*?:(.*?:.*?) ')
            barcode_id = p.findall(item[0])[0]
            if(barcode_id == ""):
                print("Got empty barcode_id")
                item = []
                continue
            if barcode in barcodes.keys():
                # res_fd = open(output_path + str(barcode) + ".fastq", "a+")
                # res_fd.write(item[0])
                # res_fd.write(item[1])
                # res_fd.write(item[2])
                # res_fd.write(item[3])
                # res_fd.close()
                barcodes[barcode].append(barcode_id + "*" + str(count - 3))
            else:
                # res_fd = open(output_path + str(barcode) + ".fastq", "a+")
                # res_fd.write(item[0])
                # res_fd.write(item[1])
                # res_fd.write(item[2])
                # res_fd.write(item[3])
                # res_fd.close()
                barcodes[barcode] = []
                barcodes[barcode].append(barcode_id + "*" + str(count - 3))
            item = []

        count += 1

    barcodes_sorted = sorted(barcodes.items(), key=lambda d:len(d[1]), reverse = True)
    barcodes_need = barcodes_sorted[0:barcodes_need_count]

    i = 0
    correct_count = 0
    error_count = 0
    for it in barcodes_sorted:
        if(i < barcodes_need_count):
            correct_count += len(it[1])
        else:
            error_count += len(it[1])
        i += 1
    fd.close()

    print("Correct rate = " + str((correct_count / (correct_count + error_count))*100) + "% (" + str(correct_count) + "/" + str(error_count) + ")")

    # print(barcodes_need)
    return barcodes_need

# type:id表转id_type表
def id_match_type(type_match_id):
    id_match_type = {}
    for item in type_match_id:
        _type = item[0]
        for i in range(len(item[1])):
            id_line = item[1][i].split('*')
            id = id_line[0]
            line = id_line[1]
            id_match_type[id] = [_type, int(line)]
    return id_match_type

# barcodes_R1: 从文件1(R1或R2)得到的type-id表
# id_match_type_R2: 从文件2得到的type-id转换后的id-type表
# input_file: 文件1输入文件（将会从此文件读取指定行输出到结果文件）
# output_path: 输出文件路径
def split_to_type_file(barcodes_R1, id_match_type_R2, file_sample_meta_data, input_file, output_path, out_suffix, barcode_line):
    found_id_count = 0
    miss_id_count = 0
    open_files = {}

    for item in barcodes_R1:
        R1_type = item[0]
        for i in range(len(item[1])):
            R1_id_line = item[1][i].split('*')
            R1_id = R1_id_line[0]
            R1_line = int(R1_id_line[1])
            if R1_id in id_match_type_R2.keys():
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
                for it in re:
                    wr_fd.write(it)

                found_id_count += 1
            else:
                print("Warning: id = " + R1_id + " is not found")
                miss_id_count += 1
                pass

    for key in open_files:
        open_files[key].close()

    return found_id_count, miss_id_count

def load_metadata(meta_file, pass_line, barcode_len):
    line_count = 0
    type_sample = {}
    fr = open(meta_file, 'r')
    for line in fr.readlines():
        line_count += 1
        if line_count <= pass_line:
            continue
        lineArr = line.strip().split()
        type_sample[lineArr[1]] = lineArr[0]
        type_sample[lineArr[1][barcode_len:] + lineArr[1][:barcode_len]] = lineArr[0]
    return type_sample

def usage(script_name):
    # python split_libraries.py L1P1_R1.fastq L1P1_R2.fastq metadata-L1P1.txt ./result 4 12 12 8 2 _R1 _R2
    print("python " + script_name +" R1_file_path R2_file_path Meta_file_path output_path barcode_line barcode_len R1_type_count R2_type_count meta_pass_line R1_output_suffix R2_output_suffix")

if __name__ == "__main__":

    starttime = datetime.datetime.now()
    global_starttime = starttime

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

    # 提取barcode类型
    barcodes_R1_need = extract_barcode(R1_file, R1_type_count, barcode_len, barcode_line)
    barcodes_R2_need = extract_barcode(R2_file, R2_type_count, barcode_len, barcode_line)

    endtime = datetime.datetime.now()
    print("extract barcode spend " + str((endtime - starttime).seconds) + " sec")

    starttime = datetime.datetime.now()

    id_match_type_R1 = id_match_type(barcodes_R1_need)
    id_match_type_R2 = id_match_type(barcodes_R2_need)

    endtime = datetime.datetime.now()
    print("trans type-id table spend " + str((endtime - starttime).seconds) + " sec")

    starttime = datetime.datetime.now()

    delete_dire(output_path) 

    os.makedirs(output_path + "/R1_R2")
    os.makedirs(output_path + "/R2_R1")

    type_sample = load_metadata(metadata_file, meta_pass_line, barcode_len)

    R1_found_id_count, R1_miss_id_count = split_to_type_file(barcodes_R1_need, id_match_type_R2, type_sample, R1_file, output_path + "/R1_R2/", R1_output_suffix, barcode_line)
    R2_found_id_count, R2_miss_id_count = split_to_type_file(barcodes_R2_need, id_match_type_R1, type_sample, R2_file, output_path + "/R2_R1/", R2_output_suffix, barcode_line)

    endtime = datetime.datetime.now()
    print("slpit to type file spend " + str((endtime - starttime).seconds) + " sec")
    print("Total time = " + str((endtime - global_starttime).seconds) + " sec")
    print("Got R1 = " + str(R1_found_id_count) + " found | " + str(R1_miss_id_count) + " miss and R2 = " + str(R2_found_id_count) + " found | " + str(R2_miss_id_count) + " miss")
