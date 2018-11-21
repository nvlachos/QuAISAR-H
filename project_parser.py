import sys
import glob
import math

def do_AR(input_AR, input_plas, output_file):
	all_ARs_in_file=[]
	samples=[]
	AR_file=open(input_AR,'r')
	line = AR_file.readline().strip()
	counter=0
	while line != '':
		#print(counter, line)
		line_sections=line.split("	")
		ar_list=line_sections[4].split(",")
		ar_dict={}
		for ar_gene in ar_list:
			gene_name=ar_gene.split(")")[0]+")"
			gene_stats=ar_gene.split(")")[1]
			ar_dict[gene_name]=gene_stats
			if gene_name not in all_ARs_in_file:
				all_ARs_in_file.append(gene_name)
				#print("Adding", gene)
		#print("1:",line_sections[0])
		#print("0:", line_sections[0], "1:", line_sections[1],"2:" , line_sections[2], "3:", line_sections[3])
		samples.append([line_sections[0], line_sections[1], line_sections[2], line_sections[3], ar_dict])
		#print("Total AR genes in sample set:", len(all_ARs_in_file)-1)
		line = AR_file.readline().strip()
	all_ARs_in_file.sort()
	if len(all_ARs_in_file) == 0:
		print("\n")
		print("Total AR genes in sample set: 0")
	else:
		print("Total AR genes in sample set:",len(all_ARs_in_file))
		#print(*all_ARs_in_file, sep = "\n")
		for gene in all_ARs_in_file:
			if gene != "No other AR genes":
				print (gene)
	print()
	AR_file.close


	#Parse plasmid summary file
	all_plasmids_in_file=[]
	plas_file=open(input_plas, 'r')
	line = plas_file.readline().strip()
	sample_p_plasmids_dict={}
	sample_f_plasmids_dict={}
	current_id=""
	counter=0
	while line != '':
		#print(counter, line)
		line_sections=line.split("	")
		#print("Current id:", current_id, ":", line_sections[0])
		if current_id == "":
			current_id = line_sections[0]+"/"+line_sections[1]
		if line_sections[0]+"/"+line_sections[1] != current_id:
			#print("New name!")
			for sample_index in range(0,len(samples)):
				#print("Looking for", current_id, ", found", samples[sample_index][0].strip())
				if current_id == samples[sample_index][0]+"/"+samples[sample_index][1].strip():
					#print(samples[sample_index], "adding", sample_f_plasmids, "and", sample_p_plasmids)
					samples[sample_index].append(sample_f_plasmids_dict)
					samples[sample_index].append(sample_p_plasmids_dict)
					break
				#else:
					#print("Sample", current_id, "does not exist")
			current_id=line_sections[0]+"/"+line_sections[1].strip()
			sample_f_plasmids_dict={}
			sample_p_plasmids_dict={}
		source_assembly=line_sections[2]
		print("Test:"+line_sections[4]+":")
		plas_perc_id=math.floor(line_sections[4])
		plas_perc_length=math.floor(100*line_sections[5].split("/")[0]/line_sections[5].split("/")[1])
		#plas_match_info="["+plas_perc_id+"/"+plas_percpercent_length+"]"
		if source_assembly == "full_assembly":
			#print("Adding:", line_sections[3], "to sample_f_plasmids")
			sample_f_plasmids_dict[line_sections[3]]="["+plas_perc_id+"/"+plas_perc_length+"]"
		elif source_assembly == "plasmid_assembly":
			#print("Adding:", line_sections[3], "to sample_p_plasmids")
			sample_p_plasmids_dict[line_sections[3]]="["+plas_perc_id+"/"+plas_perc_length+"]"
		if len(line_sections) > 1:
			if line_sections[3] not in all_plasmids_in_file:
				all_plasmids_in_file.append(line_sections[3])
				#print("Adding to project list", line_sections[3])
		#else:
			#print("Line length:", len(line_sections))
		#print()
		line = plas_file.readline().strip()
		counter=counter+1
	for sample_index in range(0,len(samples)):
		#print("Looking for", current_id, ", found", samples[sample_index][0].strip())
		if current_id == samples[sample_index][0]+"/"+samples[sample_index][1].strip():
			#print(samples[sample_index], "adding", sample_f_plasmids, "and", sample_p_plasmids)
			samples[sample_index].append(sample_f_plasmids_dict)
			samples[sample_index].append(sample_p_plasmids_dict)
			break
		#else:
			#print("Sample", current_id, "does not exist")
	all_plasmids_in_file.sort()
	if len(all_plasmids_in_file) == 1:
		#print("\n")
		print("Total plasmid replicons in sample set: 0")
	else:
		print("Total plasmid replicons in sample set:", len(all_plasmids_in_file)-1)
		print(*all_plasmids_in_file, sep= "\n")
	print()
	plas_file.close
	all_ar_and_plasmids=all_ARs_in_file+["|"]+all_plasmids_in_file
	#all_AR_to_write=all_ARs_in_file
	#all_AR_to_write.insert(0,",")
	#all_AR_to_write.insert(0,",")
	#all_AR_to_write=','.join(map(str, all_AR_to_write))
	header="id, Project__autocolour, Species__autocolour, MLST__autocolour,"
	for thing in all_ar_and_plasmids:
		header = header + " " + thing + "__autocolour,"
	header = header[:-1]
	summary_out=open(output_file, 'w')
	summary_out.write(header+'\n')
	#all_AR_to_write=all_AR_to_write[2:]
	#print("List:", all_ar_and_plasmids)
	#for sample in samples:
	#	print ("2:",sample[0])
	#return
	for sample in samples:
		sample_details=[sample[1], sample[0], sample[2], sample[3]]
		#print("pre:",sample)
		for gene in all_ar_and_plasmids:
			status=" "
			if gene == "|":
				sample_details.append(gene)
				continue
			if sample[4].get(gene):
				status=sample[4].get(gene)
			elif sample[5].get(gene):
				if sample[6].get(gene):
					status="F:"+sample[5].get(gene)+"/P:"+sample[6].get(gene)
				else:
					status="F:"+sample[5].get(gene)
			elif sample[6].get(gene):
				status="P:"+sample[6].get(gene)
			sample_details.append(status)
		#print("Post Sample check", sample_details)
		sample_details=','.join(map(str, sample_details))
		summary_out.write(sample_details+"\n")
	summary_out.close
do_AR(sys.argv[1], sys.argv[2], sys.argv[3])
