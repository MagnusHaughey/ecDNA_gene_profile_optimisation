



import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 
import seaborn as sns 
import sys
import datetime
import glob 
import scipy
import string
from pathlib import Path





def reformat(input):

	output = []
	for compartment in input.split(";"):
		output += list(np.full(int(compartment.split(",")[-1].replace("\n" , "")) , compartment.split(",")[0].replace("\n" , "")))
	return output



# Define all sequences of parameter A that were simulated, to make sure that only related data are plotted
all_A_sequences = [[20 , 30 , 40 , 60 , 80 , 100]]


for fitness_function_sequence in all_A_sequences:


	# First, find list of all parameter combinations from ./RESULTS/ directory 
	fitness_function_string = "./RESULTS/*Asequence="
	for val in fitness_function_sequence:
		fitness_function_string += "{}_".format(val)
	fitness_function_string += "*/"

	print("Searching for " + fitness_function_string)

	allDataPaths = sorted(list(set([word.split("_seed")[0] for word in glob.glob(fitness_function_string)])))


	# Plot data for each parameter combination
	for param_combo in allDataPaths:

		per_colony_data = []
		copyNumbers = []

		print("-> Reading data for parameter combination: " + param_combo)

		# List of all files at particular parameterisation and time point 
		all_files = glob.glob(param_combo + "*/*tissue.csv")

		# Option to randomly shuffle files 
		np.random.shuffle(all_files)

		# Find parameter values from directory path 
		param_combo_split = param_combo.split("_")

		k_detected = [int(word.split("=")[-1]) for word in param_combo_split if ("k=" in word)][0]
		s_detected = [float(word.split("=")[-1]) for word in param_combo_split if ("s=" in word)][0]
		C_detected = [float(word.split("=")[-1]) for word in param_combo_split if ("C=" in word)][0]
		p_detected = [float(word.split("=")[-1]) for word in param_combo_split if ("p=" in word)][0]
		q_detected = [float(word.split("=")[-1][:-1]) if (word[-1] == "/") else float(word.split("=")[-1]) for word in param_combo_split if ("q=" in word)][0]

		b_sequence_detected = param_combo.split("Bsequence=")[-1]
		b_sequence_detected = b_sequence_detected.split("_seed")[0]
		b_sequence_detected = b_sequence_detected.split("_")


		for file in all_files:

			print("\tReading data for: " + file)

			# Get patient seed from file path
			seed = [int(word.split("=")[-1]) for word in file.split('/') if "seed" in word][0]

			# Get A and B from file path 
			A_detected = [int(word.split("=")[-1]) for word in file.split("_") if ("A=" in word)][0]
			B_detected = [int(word.split("=")[-1]) for word in file.split("_") if ("B=" in word)][0]


			# List for pseudo-bulked cells
			pseudo_bulk_sample = []

			# Read file and add data to cumulative vector
			with open(file, "r") as f:
				data = f.readlines()


				# Sub-sample e.g. 100 cells and count ecDNA copy number in each 
				id_pass_filter = np.random.choice(np.arange(len(data)) , size = np.min([10 , len(data)]) , replace = False)
				for i in id_pass_filter:
					copyNumbers.append([A_detected , len(reformat(data[i]))])


				if (len(data) == 0):
					continue


				# Sub-sample small number (e.g. 1) of cells to match combiFish experiment
				subSample_size = 1
				id_pass_filter = np.random.choice(np.arange(len(data)) , size = np.min([subSample_size , len(data)]) , replace = False)


				for i in id_pass_filter:
					cell_ecDNA = reformat(data[i])

					# Compute percent of colocalization of gene A (oncogene) and gene B (passenger)
					colocalization = 100.0 * len([ec for ec in cell_ecDNA if ("A" in ec) and ("B" in ec)]) / len(cell_ecDNA)

					# Compute bulk copy numbers of gene A and gene B
					gene_A_CopyNumber = np.sum([ec.count("A") for ec in cell_ecDNA])
					gene_B_CopyNumber = np.sum([ec.count("B") for ec in cell_ecDNA])


					per_colony_data.append([k_detected , s_detected , A_detected , B_detected , C_detected , p_detected , q_detected , seed , colocalization , gene_A_CopyNumber , gene_B_CopyNumber])







		# Wrap all data into a dataframe
		per_colony_data_df = pd.DataFrame(per_colony_data , columns = ['k' , 's' , 'A' , 'B' , 'C' , 'p' , 'q' , 'seed' , 'Colocalization' , 'Gene A copy number' , 'Gene B copy number'])
		copyNumbers_df = pd.DataFrame(copyNumbers , columns = ['A' , 'Copy number'])



		# Plot colocalization of genes A and B
		plt.clf()
		plt.close()
		plt.subplots(figsize = (1+(0.5*len(fitness_function_sequence)),4.4))
		sns.boxplot(data = per_colony_data_df ,
				x = "A" ,
				y = "Colocalization" ,
				color = "#57ba96" ,
				order = fitness_function_sequence ,
				showfliers = False)

		sns.stripplot(data = per_colony_data_df ,
				x = "A" ,
				y = "Colocalization" ,
				color = "#0d6b49" ,
				edgecolor = "#4f4f4f" ,
				order = fitness_function_sequence ,
				linewidth = 0.75   ,
				size = 5 )


		plt.ylim([-5 , 105])

		plt.ylabel("Colocalization (%)" , fontsize = 15) 

		plt.tick_params(labelsize = 12 , width = 2.2)
		plt.tick_params(which = 'minor' , width = 2.2)

		plt.gca().spines[['right', 'top']].set_visible(False)
		for axis in ['bottom','left']:
			plt.gca().spines[axis].set_linewidth(2.2)

		plt.tight_layout()


		param_list = [word for word in param_combo.split("/") if ("Nmax" in word)][0]

		Path("PLOTS").mkdir(parents = True , exist_ok = True)
		outFile = param_combo.replace("RESULTS" , "PLOTS") + "_colocalization.png"

		plt.savefig(outFile , dpi = 600 , transparent = False , format = 'png')
		
		print("Plotted: " + outFile)






		# Plot ecDNA copy number distributions as a function of drug concentration
		plt.clf()
		plt.close()
		plt.subplots(figsize = (1+(0.5*len(fitness_function_sequence)),4.4))
		sns.boxplot(data = copyNumbers_df ,
				x = "A" ,
				y = "Copy number" ,
				color = "#b536b3" ,
				order = fitness_function_sequence ,
				showfliers = False)

		sns.stripplot(data = copyNumbers_df ,
				x = "A" ,
				y = "Copy number" ,
				color = "#780276" ,
				edgecolor = "#4f4f4f" ,
				order = fitness_function_sequence ,
				linewidth = 0.75   ,
				size = 5 )



		plt.ylabel("Mean copy number" , fontsize = 15) 

		plt.tick_params(labelsize = 12 , width = 2.2)
		plt.tick_params(which = 'minor' , width = 2.2)

		plt.gca().spines[['right', 'top']].set_visible(False)
		for axis in ['bottom','left']:
			plt.gca().spines[axis].set_linewidth(2.2)

		plt.tight_layout()


		param_list = [word for word in param_combo.split("/") if ("Nmax" in word)][0]

		Path("PLOTS").mkdir(parents = True , exist_ok = True)
		outFile = param_combo.replace("RESULTS" , "PLOTS") + "_meanCopyNumber.png"

		plt.savefig(outFile , dpi = 600 , transparent = False , format = 'png')
		
		print("Plotted: " + outFile)






		# Copy numbers of genes A and B
		plt.clf()
		plt.close()


		# Create mosaic and axes map for this sequence of A values (i.e. fitness functions)
		axes_map = {}
		for i in range(len(fitness_function_sequence)):
			axes_map[sorted(fitness_function_sequence)[i]] = string.ascii_uppercase[i]

		mosaic = "".join(list(axes_map.values()))
		fig , axes = plt.subplot_mosaic(mosaic , sharey = True , figsize = (len(fitness_function_sequence)*3 , 4))


		# # Construct list of mean and std for gene A and B copy numbers at each A value (i.e. each fitness curve)
		# # List is of the format [A , (mean_A_copyNumber , mean_B_copyNumber) , (std_A_copyNumber , std_B_copyNumber)]
		# gene_values_toPlot = []
		# for a_value in fitness_function_sequence:
		# 	all_gene_A_copyNumbers = per_colony_data_df[per_colony_data_df['A'] == a_value]["Gene A copy number"].to_list()
		# 	all_gene_B_copyNumbers = per_colony_data_df[per_colony_data_df['A'] == a_value]["Gene B copy number"].to_list()

		# 	mean_A = np.mean(all_gene_A_copyNumbers) if len(all_gene_A_copyNumbers) > 0 else 0
		# 	mean_B = np.mean(all_gene_A_copyNumbers) if len(all_gene_B_copyNumbers) > 0 else 0
		# 	std_A = np.std(all_gene_A_copyNumbers) if len(all_gene_A_copyNumbers) > 0 else 0
		# 	std_B = np.std(all_gene_A_copyNumbers) if len(all_gene_B_copyNumbers) > 0 else 0

		# 	gene_values_toPlot.append([a_value , (mean_A , mean_B) , (std_A , std_B)])



		# Plot 
		toPlot_x = [0.5 , 1.0 , 1.0 , 2.0 , 2.0 , 3.0 , 3.0 , 3.5]
		toPlot_y = []
		for i in range(len(fitness_function_sequence)):

			a_value = sorted(fitness_function_sequence)[i]
			b_value = b_sequence_detected[i]

			all_gene_A_copyNumbers = per_colony_data_df[per_colony_data_df['A'] == a_value]["Gene A copy number"].to_list()
			all_gene_B_copyNumbers = per_colony_data_df[per_colony_data_df['A'] == a_value]["Gene B copy number"].to_list()

			mean_A = np.mean(all_gene_A_copyNumbers) if len(all_gene_A_copyNumbers) > 0 else 0
			mean_B = np.mean(all_gene_B_copyNumbers) if len(all_gene_B_copyNumbers) > 0 else 0
			std_A = np.std(all_gene_A_copyNumbers) if len(all_gene_A_copyNumbers) > 0 else 0
			std_B = np.std(all_gene_B_copyNumbers) if len(all_gene_B_copyNumbers) > 0 else 0

			# Construct representation of ecDNA sashimi plot 
			toPlot_y = [0.0 , 0.0 , mean_A , mean_A , mean_B , mean_B , 0.0 , 0.0]


			# Plot
			axes[axes_map[a_value]].plot(toPlot_x , toPlot_y , lw = 3 , color = "black")

			
			# Highlight gene A and B region in colour
			if (a_value == sorted(fitness_function_sequence)[0]):
				axes[axes_map[a_value]].plot([1.25 , 1.75] , [mean_A , mean_A] , lw = 6 , color = "#f24324" , label = "Oncogene")
				axes[axes_map[a_value]].plot([2.25 , 2.75] , [mean_B , mean_B] , lw = 6 , color = "#6624f2" , label = "Passenger")
			else:
				axes[axes_map[a_value]].plot([1.25 , 1.75] , [mean_A , mean_A] , lw = 6 , color = "#f24324")
				axes[axes_map[a_value]].plot([2.25 , 2.75] , [mean_B , mean_B] , lw = 6 , color = "#6624f2")


			# Plot error bars on mean gene copy numbers
			A_errs = [[std_A if (mean_A - std_A) >= 0 else mean_A] , [std_A]]
			B_errs = [[std_B if (mean_B - std_B) >= 0 else mean_B] , [std_B]]
 
			axes[axes_map[a_value]].errorbar(1.5 , mean_A , yerr = A_errs , xerr = None , ecolor = "#a61900" , elinewidth = 2 , capsize = 2 , capthick = 2 , zorder = -10)
			axes[axes_map[a_value]].errorbar(2.5 , mean_B , yerr = B_errs , xerr = None , ecolor = "#3c04b5" , elinewidth = 2 , capsize = 2 , capthick = 2 , zorder = -10)


			# Formatting
			axes[axes_map[a_value]].spines[['right', 'top']].set_visible(False)
			for axis in ['bottom','left']:
				axes[axes_map[a_value]].spines[axis].set_linewidth(2.2)

			axes[axes_map[a_value]].set_xlabel(None)
			
			axes[axes_map[a_value]].tick_params(axis = 'x' , which = 'both' , bottom = False , labelbottom = False)
			axes[axes_map[a_value]].tick_params(axis = 'y' , which = 'both' , width = 2.2)

			axes[axes_map[a_value]].set_title("A={}, B={}".format(a_value, b_value) , pad = 15)



		# Legend
		handles, labels = axes["A"].get_legend_handles_labels()
		axes["A"].legend(handles, labels, loc = 'upper right' , borderpad = 1)



		axes["A"].set_ylabel("Copy number" , fontsize = 18 , labelpad = 15) 
		axes["A"].tick_params(axis = 'y' , which = 'both' , width = 2.2 , labelsize = 15)



		# labels = plt.gca().get_xticklabels()
		# plt.gca().set_xticks([word.get_position()[0] for word in labels])
		# plt.gca().set_xticklabels([word.get_text().replace(" copy number" , "") for word in labels])

		plt.tight_layout()
		plt.subplots_adjust(wspace = 0.5)

		Path("PLOTS").mkdir(parents = True , exist_ok = True)
		outFile = param_combo.replace("RESULTS" , "PLOTS") + "_copyNumbers.png"

		plt.savefig(outFile , dpi = 600 , transparent = False , format = 'png')
		
		print("Plotted: " + outFile + "\n")







