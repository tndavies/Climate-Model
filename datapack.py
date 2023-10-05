# Datapack: takes path to a datapack file (.dp) and parses it,
# exposing the comment (if present), and datapoints.
class Datapack:
	def __init__(self, datapack_name):
		self._src = datapack_name
		self.Xdata = [] 
		self.Ydata = [] 
		self.desc = ""

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
			
			x = float( arr[0] )
			y = float( arr[1].strip("\n") ) # remove line feed.
			
			self.Xdata.append(x)
			self.Ydata.append(y)

		datapack_file.close()