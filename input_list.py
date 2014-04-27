#!/usr/bin/python

# function to get the column number of a variable
def parse_input(data_file, variable):
	table_dir = "/net/edwards/ta263/agetables/"
	file_path = table_dir + data_file
	num_entries = 25 # 25 for just line intensities, 16 for line ratios
	test_bad_row = 12 # if there is a 'NaN' in this row, then we don't pull data from any columns in this row. 12 for intensities, 14 for ratios. 
	bad_data_value = '-100.0'

# get the first line of the file so we can see which column # 'variable' refers to
	f = open(file_path, 'r')
	output_data = []
	header_line = f.readline()

# figure out which column 'variable' belongs to
# loop through the words in this sentence to do this
	i = 1
	col_num = 'NaN'
	for word in header_line.split():
		if word==variable: 
			col_num = i
			break
		else:
			i = i+1

# THIRD loop through the rest of the file and write 
	for line in f:
		this_entries = len(line.split())
		# loop through lines with data in them, copy values from good fibres (definition of good fibre = number in col 13)
		if (this_entries == num_entries) and not 'NaN' in line.split()[test_bad_row]:
			output_data.append(float(line.split()[col_num-1]))
		# else, if there is a bad value in that fibre, write -100 to the file
		elif (this_entries == num_entries) and 'NaN' in line.split()[test_bad_row]:
			output_data.append(float(bad_data_value))
		# else, move on

	f.close()
	
	return output_data

