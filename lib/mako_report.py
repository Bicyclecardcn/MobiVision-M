import json
from mako.template import Template
import base64
import re 
import numpy as np
import os

def shorten_text(i_txt:str, shorten_len:int):
    if len(i_txt) <= shorten_len:
        return i_txt
    else:
        return i_txt[:shorten_len] + "..."

def show_number(n, to_int=False, to_sci=False):
    if n > 1e6 and to_sci:
        return_n = "{:.2e}".format(n)
    elif n % 1 == 0 or to_int:
        return_n = "{:,.0f}".format(n)
    else:
        return_n = "{:,.2f}".format(n)
    return return_n

def get_multi_select(json_data:dict):
    return_str = ""
    n = 1
    for i in json_data["species_info"].keys():
        return_str += '<option value="%s">%s</option>\n' %(show_number(n), i)
        i += 1
    return return_str

def get_input_info(json_data:dict):
    return_str = ""
    export_dict = {"raw_reads":"Number of Raw Reads", 
                   "valid_reads":"Number of Valid Reads", 
                   "fraction_of_host_reads":"Fraction of Host Reads", 
                   "Fraction_of_reads":"Fraction of Valid Reads", 
                   "saturation_rate":"Sequencing Saturation", 
                   "Q30_in_valid":"Q30 Bases in Barcode+UMI", 
                   "Q30_in_all":"Q30 Bases in RNA Read"}
    for i in json_data["input"].keys():
        if i in export_dict.keys():
            show_data = json_data["input"][i]
            if i in ["Fraction_of_reads", "saturation_rate", "fraction_of_host_reads", "Q30_in_valid", "Q30_in_all"]:
                show_data = show_number(show_data*100) + "%"
            else:
                show_data = show_number(show_data)
            return_str += '<div class="data"><span class="label" en="%s" zh=""></span>\n' %(export_dict[i])+ \
                        '<div class="data-num">%s</div></div>\n' %(show_data)
    return return_str

def get_multi_species1(json_data:dict):
    ###cell number
    return_str = '<details id="cell_number_details">'
    return_str += '<summary class="data"><span class="label sub-data-label" en="Estimated Number of Microbes" zh=""></span>' + \
                '<div class="data-num">%s</div></summary>' %(show_number(json_data["species_info"]["all"]["cell_number"], to_int=True))
    for i in json_data["species_info"].keys():
        if i != "all" and i != "unknown":
            #if i == "unknown":
            #    show_tag = "Number of barcodes with > 1 Microbes"
            #else:
            #    show_tag = "Estimated Number of Microbes (%s)" %(i)
            show_data = show_number(json_data["species_info"][i]["cell_number"])
            return_str += '<div class="subdata"><span class="label" en="%s" zh=""></span>' %(i) + \
                        '<div class="data-num">%s</div></div>' %(show_data)
    for i in json_data["species_info"].keys():
        if i == "unknown":
            #show_tag = "Number of barcodes with > 1 Microbes"
            show_tag = "Multiplet "
            show_data = show_number(json_data["species_info"][i]["cell_number"])
            return_str += '<div class="subdata"><span class="label" en="%s" zh=""></span>' %(show_tag) + \
                        '<div class="data-num">%s</div></div>' %(show_data)
    return_str += '</details>'
    ###mean reads
    return_str += '<details id="mean_reads_details">'
    show_data = show_number(json_data["species_info"]["all"]["mean_reads_count"], to_int=True)
    return_str += '<summary class="data"><span class="label sub-data-label" en="Mean reads per Microbe" zh=""></span>' + \
                '<div class="data-num">%s</div></summary>' %(show_data)
    for i in json_data["species_info"].keys():
        #if i != "unknown" and i != "all":
        if i != "all":
            show_data = json_data["species_info"][i]["mean_reads_count"]
            if isinstance(show_data, float):
                #show_data = "%.2f" %(show_data)
                show_data = show_number(show_data, to_int=True)
            if i == "unknown":
                show_tag = "Multiplet "
            else:
                show_tag = i
            return_str += '<div class="subdata"><span class="label" en="%s" zh=""></span>' %(show_tag) + \
                        '<div class="data-num">%s</div></div>' %(show_data)
    return_str += "</details>"
    ###Median gene
    return_str += "<details id='median_gene_details'>"
    show_data = show_number(json_data["species_info"]["all"]["median_gene_count"], to_int=True)
    return_str += '<summary class="data"><span class="label sub-data-label" en="Median Genes per Microbe" zh=""></span>' + \
                '<div class="data-num">%s</div></summary>' %(show_data)
    for i in json_data["species_info"].keys():
        #if i != "unknown" and i != "all":
        if i != "all":
            show_data = show_number(json_data["species_info"][i]["median_gene_count"])
            if i == "unknown":
                show_tag = "Multiplet "
            else:
                show_tag = i
            return_str += '<div class="subdata"><span class="label" en="%s" zh=""></span>' %(show_tag) + \
                        '<div class="data-num">%s</div></div>' %(show_data)
    return_str += "</details>"
    ###median UMI
    return_str += "<details id='median_UMI_details'>"
    show_data = show_number(json_data["species_info"]["all"]["median_UMI_count"], to_int=True)
    return_str += '<summary class="data"><span class="label sub-data-label" en="Median UMI per Microbe" zh=""></span>' + \
                '<div class="data-num">%s</div></summary>' %(show_data)
    for i in json_data["species_info"].keys():
        #if i != "unknown" and i != "all":
        if i != "all":
            show_data = show_number(json_data["species_info"][i]["median_UMI_count"])
            if i == "unknown":
                show_tag = "Multiplet "
            else:
                show_tag = i
            return_str += '<div class="subdata"><span class="label" en="%s" zh=""></span>' %(show_tag) + \
                        '<div class="data-num">%s</div></div>' %(show_data)
    return_str += "</details>"
    ###gene detected
    return_str += "<details id='gene_detected_details'>"
    show_data = show_number(json_data["species_info"]["all"]["gene_detected"])
    return_str += '<summary class="data"><span class="label sub-data-label" en="Total Gene Detected" zh=""></span>' + \
                '<div class="data-num">%s</div></summary>' %(show_data)
    for i in json_data["species_info"].keys():
        #if i != "unknown" and i != "all":
        if i != "all" and i != "unknown":
            show_data = show_number(json_data["species_info"][i]["gene_detected"])
            return_str += '<div class="subdata"><span class="label" en="%s" zh=""></span>' %(i) + \
                        '<div class="data-num">%s</div></div>' %(show_data)
    return_str += "</details>"
    return return_str

def get_multi_species2(json_data:dict):
    return_str = "<details>"
    return_str += '<summary class="data"><span class="label sub-data-label" en="Mapping Information" zh=""></span>' + \
                '<div class="data-num"></div></summary>'
    export_dict = {"map_rate":"Reads Mapped to Genome", "confidently_map_rate":"Reads Mapped Confidently to Genome"}
    for i in json_data["species_info"].keys():
        if i != "unknown":
            for j in json_data["species_info"][i].keys():
                if j in export_dict.keys():
                    show_data = json_data["species_info"][i][j]
                    if isinstance(show_data, float):
                        show_data = show_number(show_data*100)
                    return_str += '<div class="subdata"><span class="label" en="%s (%s)" zh="">' %(export_dict[j], i) + \
                                '</span><div class="data-num">%s</div></div>\n' %(show_data + "%")
    return_str += "</details>"
    return return_str

def get_js(js_file:str):
    return_str = ""
    with open(js_file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            return_str += line
    return return_str

def image_to_base64(image_dir:str):
    base64_data = ""
    if image_dir is not None:
        if os.path.exists(image_dir):
            with open(image_dir, "rb") as f:
                base64_data = str(base64.b64encode(f.read()))
                if base64_data[:2] == "b'":
                    base64_data = "'data:image/png;base64," + base64_data[2:]
    return base64_data

def find_first_difference(str1, str2):
    min_length = min(len(str1), len(str2))
    for i in range(min_length):
        if str1[i] != str2[i]:
            echo = min_length - i - 1
            if echo > 3:
                echo = 3
            return i, str1[i+echo], str2[i+echo]
    if len(str1) > len(str2):
        return min_length, str1[min_length:(len(str1) - min_length)], ""
    elif len(str1) < len(str2):
        return min_length, "", str2[min_length:(len(str2) - min_length)]
    return -1, "", ""

def shorten_species_names(json_data:dict):
    name_dict = {}
    for i in json_data.keys():
        if len(i) > 20:
            tmp_name = i[:20]
        else:
            tmp_name = i
        #while tmp_name in name_dict.keys():
        #    prefix = tmp_name
        #    start, df1, df2 = find_first_difference(i, name_dict[tmp_name])
        #    if start == -1:
        #        print("Error! Species names the same: %s and %s." %(i, name_dict[tmp_name]))
        #        sys.exit()
        #    else:
        #        if start <= 13:
        #            update_name1 = prefix + df2
        #            name_dict[update_name1] = name_dict[tmp_name]
        #            del name_dict[tmp_name]
        #            tmp_name = prefix + df1
        #        else:
        #            update_name1 = prefix + "..." + df2
        #            name_dict[update_name1] = name_dict[tmp_name]
        #            del name_dict[tmp_name]
        #            tmp_name = prefix + "..." + df1
        name_dict[tmp_name] = i
    return_dict = {}
    for i in name_dict.keys():
        return_dict[i] = json_data[name_dict[i]]
    return return_dict

class ExportReport():
    def __init__(self, template_file:str, json_file:str, output_file:str, jquery:str, plotly:str, favicon_file:str, web_logo:str, web_back:str) -> None:
        self.template_file = template_file
        self.json_file = json_file
        self.output_file = output_file
        self.jquery_script_file = jquery
        self.plotly_script_file = plotly
        self.favicon_file = favicon_file
        self.web_logo_file = web_logo
        self.web_back_file = web_back

    def process(self):
        with open(self.json_file, "r") as f:
            json_data = json.load(f)
        json_data["species_info"] = shorten_species_names(json_data["species_info"])
        web_keywords = "墨卓生物科技,数字PCR,单细胞测序平台,微流体检测芯片"
        web_description = "墨卓生物是中国首家多平台新型分子诊断公司。基于哈佛大学最前沿的第三代微流控技术专利矩阵，墨卓生物致力于打造更精准、更快速、更实用的新型分子诊断与生命科学研究平台。升级诊断与检测质量，提高生命科学研究水平，为科研人员与医生、患者带来更高效的解决方案。"
        ##json_data["multi_select"] = "" ###https://blog.csdn.net/nayi_224/article/details/86136785, https://github.com/wenzhixin/multiple-select
        web_viewport = "width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no"
        web_favicon = image_to_base64(self.favicon_file)
        web_log = image_to_base64(self.web_logo_file)
        web_back = image_to_base64(self.web_back_file)
        jquery_script = get_js(self.jquery_script_file)
        plotly_script = get_js(self.plotly_script_file)
        #multi_select_info = get_multi_select(json_data)
        multi_species_info1 = get_multi_species1(json_data)
        input_info = get_input_info(json_data)
        multi_species_info2 = get_multi_species2(json_data)

        mytemplate = Template(filename=self.template_file)
        tmp_foot = "MobiVision-M %s finished at %s. <br>" %(json_data["sample"]["version_detail"], json_data["sample"]["run_time"])
        tmp_args = "With specific arguments:"
        run_pattern = {"-t":r"-t (\d+)", 
                       "-s":r"-s (\S+)", 
                       "-c":r"-c (\S+)",
                       "--call_mtx":r"--call_mtx (\S+)",
                       "--sample_ID":r"--sample_ID (\S+)",
                       "--cellnumber":r"--cellnumber (\d+)", 
                       "--UMI_adjust":r"--UMI_adjust (\S+)", 
                       "--multiplet_method":r"--multiplet_method (\S+)", 
                       "--host_remove":r"--host_remove (\S+)", 
                       "--host_reference":r"--host_reference (\S+)",
                       "--config":r"--config (\S+)", 
                       "--hard_filter":r"--hard_fileter (\S+)",
                       "--kit":r"--kit (\S+)"}
        switch_arg = [" --cr2.2 ", " --nosecondary ", " --keep_bam ", " --keep_unmap_reads ", " --qc_only ", " --test_run "]
        for i in switch_arg:
            if i in json_data["sample"]["run_cmd"]:
                tmp_args += i
        for i in run_pattern.keys():
            match = re.search(run_pattern[i], json_data["sample"]["run_cmd"])
            if match:
                tmp_args += " %s %s" %(i, match.group(1))
        foot_height = str((np.ceil(len(tmp_args) / 150) + 1) * 20) + "px"
        if json_data["all_plot"]["UMAP_count"] != {}:
            UMAP_display = "block"
        else:
            UMAP_display = "none"
        with open(self.output_file, "w") as f:
            f.write(mytemplate.render(input_json=json_data["all_plot"], 
                                    Sample_ID=shorten_text(i_txt=json_data["sample"]["id"], shorten_len=25),
                                    reference=shorten_text(i_txt=json_data["sample"]["reference"], shorten_len=25), 
                                    kit=shorten_text(json_data["sample"]["kit"], shorten_len=25), 
                                    version=shorten_text(json_data["sample"]["version"], shorten_len=25), 
                                    all_cell_number=show_number(json_data["species_info"]["all"]["cell_number"], to_sci=True),
                                    all_median_gene=show_number(json_data["species_info"]["all"]["median_gene_count"], to_sci=True), 
                                    all_mean_reads=show_number(json_data["species_info"]["all"]["mean_reads_count"], to_sci=True), 
                                    all_median_UMI=show_number(json_data["species_info"]["all"]["median_UMI_count"], to_sci=True),  
                                    input_info=input_info,
                                    multi_species_info1=multi_species_info1,
                                    multi_species_info2=multi_species_info2,
                                    web_keywords=web_keywords, 
                                    web_description=web_description, 
                                    web_viewport=web_viewport, 
                                    web_favicon=web_favicon, 
                                    web_logo=web_log,
                                    web_back=web_back, 
                                    jquery_js=jquery_script, 
                                    plotly_js=plotly_script, 
                                    web_title=json_data["sample"]["id"], 
                                    foot=tmp_foot,
                                    args=tmp_args, 
                                    foot_height=foot_height, 
                                    run_cmd=json_data["sample"]["run_cmd"], UMAP_display=UMAP_display))

    
if __name__ == "__main__":
    import argparse
    parse = argparse.ArgumentParser()
    parse.add_argument('--template_file', type=str, help='The path of html template.')
    parse.add_argument('--json_file', type=str, help="The path of reference json.")
    parse.add_argument('-o', '--output_file', type=str, help="The path of output report.")
    parse.add_argument('--jquery_script_file', type=str, help="The path of jquery script file.")
    parse.add_argument('--plotly_script_file', type=str, help="The path of plotly.js script file.")
    parse.add_argument('--favicon_file', type=str, help="The path of favicon file.")
    parse.add_argument('--web_logo_file', type=str, help="The path of web logo file.")
    parse.add_argument('--web_back_file', type=str, help="The path of web back file.")
    args = parse.parse_args()
    p = ExportReport(template_file=args.template_file, 
                     json_file=args.json_file, 
                     output_file=args.output_file, 
                     jquery=args.jquery_script_file, 
                     plotly=args.plotly_script_file,
                     favicon_file=args.favicon_file, 
                     web_logo=args.web_logo_file, 
                     web_back=args.web_back_file)
    p.process()
