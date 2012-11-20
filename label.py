
class instance_label(object):

	def __init__(self):

		self.examplars = []
		self.instance_label = {}

	def read_examplars(self, filename):

		f = file(filename)
		lines = f.readline().strip().split()
		for i in xrange(len(lines)):
			self.examplars.append(lines[i])
	
		#print self.examplars
	
	def mark_instance_label(self, filename):

		f = file(filename)		
		
		for line in f:
			lines = line.strip().split()
			
			if not lines[0] in self.examplars and lines[1] in self.examplars:
				if not lines[0] in self.instance_label:
					self.instance_label[lines[0]] = [lines[1], lines[2]]
				else:
					if float(lines[2]) > float(self.instance_label[lines[0]][1]):
						self.instance_label[lines[0]][0] = lines[1]
						self.instance_label[lines[0]][1] = lines[2]

		
		for key, value in self.instance_label.items():
			print key, value
							
		

if __name__ == "__main__":
	examplars = "examplars"

	il = instance_label()
	il.read_examplars(examplars)

	filename = "output"
	il.mark_instance_label(filename)

	
