import sys
import os

import ROOT

def loop_merge(source_dir, dir_name, temp_dir, n_selections, append, final_name):
    last_failed = False
    n_events = 0
    data_tree_list_per_selection = [[] for i in range(n_selections)]
    file_list_per_selection = [[] for i in range(n_selections)]
    t_list_per_selection = [ROOT.TList() for i in range(n_selections)]
    merge_count = 0

    fail_loop = False
    finished_loop = True
    while 1:
        for i in range(n_selections):
            fname = source_dir + str(dir_name) + "/out%i.root" % i
            if not os.path.isfile(fname):
                if last_failed:
                    fail_loop = True
                    finished_loop = False
                    print "END OF DATA"
                    break
                else:
                    finished_loop = False
                    last_failed = True
                    dir_name += 1
                    break
            else:
                if last_failed:
                    print "MISSING DIRECTIORY %i " % (dir_name-1)
                last_failed = False

            file_list_per_selection[i].append(ROOT.TFile(fname))
            # check if file is corrupted
            try:
                file_list_per_selection[i][-1].Get("truth").GetEntries()
            except AttributeError:
                print "CORRUPTED FILE: %s" % fname
                if i != 0:
                    raise
                assert(i == 0)
                file_list_per_selection[i].pop()
                dir_name += 1
                finished_loop = False
                break

            data_tree_list_per_selection[i].append( file_list_per_selection[i][-1].Get("data") )
            t_list_per_selection[i].Add(data_tree_list_per_selection[i][-1])
        if not finished_loop:
            if fail_loop:
                break 
            else:
                continue

        merge_count += 1
        max_n_merge = 400 / n_selections
        if (dir_name + 1) % max((max_n_merge / 10),1) == 0:
            print "Processed %i files" % ((dir_name + 1) * n_selections)
        if merge_count % max(max_n_merge,1) == 0:
            break
        dir_name += 1

    if append:
        print "LOAD PREVIOUS AND MERGE AND WRITE" 
        for i in range(n_selections):
            temp_file_a = temp_dir + "temp%i" % i + final_name
            last_file = ROOT.TFile(temp_file_a)
            data_tree = last_file.Get("data")
            t_list_per_selection[i].Add(data_tree)

            header = file_list_per_selection[i][-1].Get("header")

            temp_file_b = temp_dir + "temp%iB" % i + final_name
            merge_and_write(t_list_per_selection[i], header, temp_file_b)

            last_file.Close()
    else:
        print "MERGE AND WRITE"
        for i in range(n_selections):
            header = file_list_per_selection[i][-1].Get("header")
            temp_file_b = temp_dir + "temp%iB" % i + final_name
            merge_and_write(t_list_per_selection[i], header, temp_file_b)

    [[f.Close() for f in lst] for lst in file_list_per_selection]
    for i in range(n_selections):
        temp_file_a = temp_dir + "temp%i" % i + final_name
        temp_file_b = temp_dir + "temp%iB" % i + final_name
        os.rename(temp_file_b, temp_file_a)

    if last_failed:
        return (0, False)
    else: 
        return (dir_name+1, True)

def merge_and_write(t_list, header, f_out_name):
    f_out = ROOT.TFile(f_out_name, "RECREATE")
    f_out.cd()
    merged = ROOT.TTree.MergeTrees(t_list)
    merged.SetName("data")

    merged.Write()
    header.Write()
    f_out.Close()

def main(source_dir, output_dir, identifier, temp_dir, n_selections):
    append = False
    dir_name = 0
    while 1:
        (dir_name, append) = loop_merge(source_dir, dir_name, temp_dir, n_selections, append, identifier)
        if not append:
            break
    for i in range(n_selections):
        os.rename(temp_dir + "temp%i" % i + final_name, output_dir + identifier + "%i.root" % i)
    print "N DATA: %i" % n_data

if __name__=="__main__":
    source_dir = sys.argv[1]
    output_dir = sys.argv[2]
    output_identifier = sys.argv[3]
    temp_dir = sys.argv[4]
    n_selections = int(sys.argv[5])
   
    main(source_dir, output_dir, output_identifier, temp_dir, n_selections) 

