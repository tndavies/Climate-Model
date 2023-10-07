# Datapack: takes path to a datapack file (.dp) and parses it,
# exposing the comment (if present), and datapoints.
class Datapack:
	def __init__(self, datapack_name):
		self._src = datapack_name
		self.desc = ""
		self.xs = [] 
		self.ys = [] 

		self.load()

	def load(self):
		datapack_file = open(self._src) # open datapack for reading.
		contents = datapack_file.readlines()

		# parse the comment, if present.
		if(contents[0][0] == '#'):
			self.desc = contents.pop(0).strip("\n") # cache & remove the comment from the list, without the line feed.

		# parse the dataset
		for s in contents:
			arr = s.split(",")
			
			try:
				x = float( arr[0] )
				y = float( arr[1].strip("\n") ) # remove line feed.
			
				self.xs.append(x)
				self.ys.append(y)
			except:
				# print("Warning: Couldn't parse a datapoint, missing it out.")
				pass

		datapack_file.close()