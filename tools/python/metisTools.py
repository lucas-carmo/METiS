import numpy as np
import matplotlib.pyplot as plt

class metis():
    def __init__(self, inputFile=None, outputFile=None):
        self.fowt = {'floater': {'morisonElements': []}}
        self.envir = {'nodes': {'ID': [], 'coord': []}}

        # Check if is string
        if inputFile is not None:
            self.readInputFile(inputFile)

        if outputFile is not None:
            self.readOutputFile(outputFile)

    def readOutputFile(self, flPath):
        if not isinstance(flPath, str):
            raise ValueError('The input must be a string with the path to a single METiS out file')

        data = np.genfromtxt(flPath, names=True)
        self.out = {name.lower(): data[name] for name in data.dtype.names}

    def readInputFile(self, flPath):
        if not isinstance(flPath, str):
            raise ValueError('The input must be a string with the path to a single METiS input file')


        with open(flPath, 'r') as file:
            lines = file.readlines()

        lines = [line.split('%', 1)[0].strip() for line in lines]  # Remove comments
        lines = [line for line in lines if line]  # Remove empty lines

        while lines:
            line = lines.pop(0)
            if line.lower().startswith('nodes'):
                while lines and not lines[0].lower().startswith('end'):
                    line = lines.pop(0)
                    data = line.split(',')
                    if len(data) != 4:
                        raise ValueError('Wrong number of parameters.')
                    self.envir['nodes']['ID'].append(float(data[0]))
                    self.envir['nodes']['coord'].append([float(data[1]), float(data[2]), float(data[3])])
            elif line.lower().startswith('floatercog'):
                data = self._getData(line).split(',')
                if len(data) != 3:
                    raise ValueError('Unable to read the CoG. Wrong number of parameters.')
                self.fowt['floater']['cog'] = [float(data[0]), float(data[1]), float(data[2])]
            elif line.lower().startswith('morison_circ'):
                while lines and not lines[0].lower().startswith('end'):
                    line = lines.pop(0)
                    data = line.split()
                    if len(data) != 11:
                        raise ValueError('Wrong number of parameters.')
                    if not self.envir['nodes']['ID']:
                        raise ValueError('Nodes should be specified before Morison Elements.')
                    morison_element = {
                        'node1': self.envir['nodes']['coord'][self.envir['nodes']['ID'].index(float(data[0]))],
                        'node2': self.envir['nodes']['coord'][self.envir['nodes']['ID'].index(float(data[1]))],
                        'diam': float(data[2]),
                        'CD': float(data[3]),
                        'CM': float(data[4]),
                        'numIntPoints': int(data[5]),
                        'botDiam': float(data[6]),
                        'topDiam': float(data[7]),
                        'axialCD': float(data[8]),
                        'axialCa': float(data[9]),
                        'botPressFlag': int(data[10])
                    }
                    self.fowt['floater']['morisonElements'].append(morison_element)
    
    def view_fowt(self, ax=None):
        fowt = self.fowt
        morison_elements = fowt['floater']['morisonElements']

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        for element in morison_elements:
            D = element['diam']
            n1 = np.array(element['node1'])
            n2 = np.array(element['node2'])

            d_theta = np.pi / 20
            theta = np.arange(0, 2*np.pi, d_theta)

            x1 = (D / 2) * np.cos(theta)
            y1 = (D / 2) * np.sin(theta)

            x2 = x1
            y2 = y1

            zvec = (n2 - n1) / np.linalg.norm(n2 - n1)
            if np.array_equal(zvec, np.array([0, 0, 1])):
                xvec = np.array([1, 0, 0])
                yvec = np.array([0, 1, 0])
            else:
                yvec = np.cross(np.array([0, 0, 1]), zvec) / np.linalg.norm(np.cross(np.array([0, 0, 1]), zvec))
                xvec = np.cross(yvec, zvec) / np.linalg.norm(np.cross(yvec, zvec))

            P1_local = np.array([x1, y1, np.zeros_like(x1)])
            P2_local = np.array([x2, y2, np.full_like(x2, np.linalg.norm(n2 - n1))])

            M_local2global = np.array([xvec, yvec, zvec]).T

            P1_global = M_local2global @ P1_local + n1[:, None]
            P2_global = M_local2global @ P2_local + n1[:, None]

            npoints = element['numIntPoints']
            X = np.zeros((npoints, len(theta)))
            Y = np.zeros_like(X)
            Z = np.zeros_like(X)
            for jj in range(len(theta)):
                X[:, jj] = np.linspace(P1_global[0, jj], P2_global[0, jj], npoints)
                Y[:, jj] = np.linspace(P1_global[1, jj], P2_global[1, jj], npoints)
                Z[:, jj] = np.linspace(P1_global[2, jj], P2_global[2, jj], npoints)

            color_matrix = np.ones_like(Z)
            color_matrix[Z > 0] = 2

            ax.plot_surface(X, Y, Z, facecolors=plt.cm.coolwarm(color_matrix), shade=False)
            ax.plot_surface(X, Y, np.zeros_like(Z), color='b', alpha=0.1)
            ax.plot_surface(X, Y, np.full_like(Z, np.linalg.norm(n2 - n1)), color='b', alpha=0.1)

            ax.scatter(*n1, color='r', s=100)
            ax.scatter(*n2, color='r', s=100)

    def _getData(self, in_str):
        # Check if input string is empty
        if not in_str:
            raise ValueError('Empty string passed to get_data()')

        str_tokenized = in_str.split()

        # If str_tokenized has only one element, i.e. only the keyword, then return an empty string
        if len(str_tokenized) == 1:
            return str_tokenized[0]

        out_str = ''
        for s in str_tokenized[1:]:
            out_str += s + '\t'

        return out_str        

# function keyword = getKeyword(inStr)
# 	% Get the first part of the string, the one before the first '\t' or white-space
#     keyword = strtok(inStr);
# end

# function outStr = getData(inStr)
# 	% Check if input string is empty
#     if isempty(inStr)
# 		error('Empty string passed to getData()');
#     end    

# 	str_tokenized = strsplit(inStr);

# 	% If str_tokenized has only one element, i.e. only the keyword, then return an empty string
#     if numel(str_tokenized) == 1
#         outStr = str_tokenized{1};
# 		return
#     end

# 	outStr = '';
#     for ii = 2 : numel(str_tokenized)		
# 		outStr = [outStr , str_tokenized{ii}, sprintf('\t')];
#     end
# end

# % Verify whether a string contains a comment, marked by a '%'
# function bool = thereIsCommentInString(inStr)
#     bool = contains(inStr, '%');
# end

# % Verify whether a string has content, i.e. if:
# % 1) it is not empty
# % 2) it is not just white spaces or tabs
# % 3) it does not start with a comment mark ('%')
# function bool = hasContent(inStr)	
#     inStr = inStr(~isspace(inStr));

#     if (isempty(inStr))
# 		bool = 0;
#         return
#     end
    
#     if inStr(1) == '%'
#         bool = 0;
#         return
#     end
    
#     bool = 1;
# end

# function outStr = removeComments(inStr)    
# 	outStr = inStr(1:strfind(inStr, '%')-1);
# end
    


