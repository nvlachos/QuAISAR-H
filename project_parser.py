import sys
import glob
import math

def do_AR(input_csstar_AR, input_plas, output_file, input_srst2_AR):
	all_ARs_in_file=[]
	samples=[]
	csstar_file=open(input_csstar_AR,'r')
	csstar_line = csstar_file.readline().strip()
	counter=0
	print("Start")
	while csstar_line != '':
		print(counter, csstar_line)
		print("Start csstar loop")
		csstar_line_sections=csstar_line.split("	")
		ar_list=csstar_line_sections[4].split(",")
		ar_dict={}
		for ar_gene in ar_list:
			gene_name=ar_gene.split("[")[0]
			if gene_name == "No AR genes discovered":
				gene_stats="[0/0]C"
			else:
				gene_stats="["+ar_gene.split("[")[1]+"C"
			ar_dict[gene_name]=gene_stats
			if gene_name not in all_ARs_in_file:
				all_ARs_in_file.append(gene_name)
				print("Adding", gene_name)
		srst2_file=open(input_srst2_AR,'r')
		srst2_line=srst2_file.readline().strip()
		counter=0
		while srst2_line != '':
			if counter == 4:
				exit()
			print("Start srst2 loop")
			spot_count=0
			#for k, v in ar_dict.items():
			#	print(spot_count)
			#	print(k, v)
			#	spot_count+=1
			print("Checking", srst2_line)
			srst2_line_sections=srst2_line.split("	")
			if csstar_line_sections[0] == srst2_line_sections[0] and csstar_line_sections[1] == srst2_line_sections[1]:
				print("Found",  csstar_line_sections[1], "in srst2 summary file")
				srst2_ar_list=srst2_line_sections[2].split(",")
				for srst2_ar_gene in srst2_ar_list:
					gene_name=srst2_ar_gene.split("[")[0]
					gene_stats="["+srst2_ar_gene.split("[")[1]+"S"
					print("Looking up", gene_name, "in csstar dic")
					if ar_dict.get(gene_name):
						if ar_dict.get(gene_name) != "No Other AR genes":
							print("Found", gene_name, "in both outputs")
							print("New value: "+ar_dict.get(gene_name)+":"+gene_stats)
							ar_dict[gene_name]=""+ar_dict.get(gene_name)+":"+gene_stats
						#else:
							#print("No AR found in csstar for", gene_name)
						#	:
					else:
						print("New gene", gene_name,"found in srst2")
						ar_dict[gene_name]=gene_stats
					if gene_name not in all_ARs_in_file:
						all_ARs_in_file.append(gene_name)
				break
			#else:
			#	#print(csstar_line_sections[1], "C does not equal S", srst2_line_sections[1])
			#	:
			srst2_line=srst2_file.readline().strip()
			counter+=1
		srst2_file.close()
		#for k, v in ar_dict.items():
		#	print(k, v)
		#print("1:",csstar_line_sections[0])
		#print("0:", csstar_line_sections[0], "1:", csstar_line_sections[1],"2:" , csstar_line_sections[2], "3:", csstar_line_sections[3])
		samples.append([csstar_line_sections[0], csstar_line_sections[1], csstar_line_sections[2], csstar_line_sections[3], ar_dict])
		#print("Total AR genes in sample set:", len(all_ARs_in_file)-1)
		csstar_line = csstar_file.readline().strip()
	csstar_file.close
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
		plasmid_line_sections=line.split("	")
		#print("Current id:", current_id, ":", plasmid_line_sections[0])
		if current_id == "":
			current_id = plasmid_line_sections[0]+"/"+plasmid_line_sections[1]
		if plasmid_line_sections[0]+"/"+plasmid_line_sections[1] != current_id:
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
			current_id=plasmid_line_sections[0]+"/"+plasmid_line_sections[1].strip()
			sample_f_plasmids_dict={}
			sample_p_plasmids_dict={}
		source_assembly=plasmid_line_sections[2]
		#print("Test:"+plasmid_line_sections[4]+":")
		if plasmid_line_sections[4].find("_contigs-") >= 0:
			line = plas_file.readline().strip()
			continue

		plas_perc_id=math.floor(float(plasmid_line_sections[4]))
		#print("testing:", plasmid_line_sections[5].split("/")[0], plasmid_line_sections[5].split("/")[1])
		plas_perc_length=(100*int(plasmid_line_sections[5].split("/")[1])//int(plasmid_line_sections[5].split("/")[0]))
		#plas_match_info="["+plas_perc_id+"/"+plas_percpercent_length+"]"
		if source_assembly == "full_assembly":
			#print("Adding:", plasmid_line_sections[3], "to sample_f_plasmids")
			sample_f_plasmids_dict[plasmid_line_sections[3]]="["+str(plas_perc_id)+"/"+str(plas_perc_length)+"]"
		elif source_assembly == "plasmid_assembly":
			#print("Adding:", plasmid_line_sections[3], "to sample_p_plasmids")
			sample_p_plasmids_dict[plasmid_line_sections[3]]="["+str(plas_perc_id)+"/"+str(plas_perc_length)+"]"
		if len(plasmid_line_sections) > 1:
			if plasmid_line_sections[3] not in all_plasmids_in_file:
				all_plasmids_in_file.append(plasmid_line_sections[3])
				#print("Adding to project list", plasmid_line_sections[3])
		#else:
			#print("Line length:", len(csstar_plasmid_line_sections))
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
	if len(all_plasmids_in_file) == 0:
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
		print("pre:",sample)
		for gene in all_ar_and_plasmids:
			status=" "
			if gene == "|":
				sample_details.append(gene)
				continue
			if sample[4].get(gene):
				status=sample[4].get(gene)
			elif sample[5].get(gene):
				if sample[6].get(gene):
					status="F:"+sample[5].get(gene)+";P:"+sample[6].get(gene)
				else:
					status="F:"+sample[5].get(gene)
			elif sample[6].get(gene):
				status="P:"+sample[6].get(gene)
			sample_details.append(status)
		#print("Post Sample check", sample_details)
		sample_details=','.join(map(str, sample_details))
		summary_out.write(sample_details+"\n")
	summary_out.close
#print (sys.argv[1:])
print("Parsing AR files ...\n")
do_AR(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
