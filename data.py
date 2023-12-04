import csv

def ImportTemps(path):
	with open(path, 'r') as file:
	    parser = csv.reader(file)
	    next(parser) # advance from header row

	    years, temps = [], []
	    for row in parser:
	    	years.append(int(row[0]))
	    	temps.append(float(row[1]))

	return years, temps

def ImportCO2(path):
	with open(path, 'r') as file:
	    parser = csv.reader(file)
	    next(parser) # advance from header row

	    years, co2_ppm = [], []
	    for row in parser:
	    	years.append(int(row[0]))
	    	co2_ppm.append(float(row[1]))

	return years, co2_ppm
