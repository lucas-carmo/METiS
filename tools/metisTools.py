#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

def readOutFl(flPath):
    mts = pd.read_csv(flPath, delim_whitespace = True)
    mts.columns = mts.columns.str.lower()
    return mts

def readInputFl(flPath):    
    fl = open(flPath, 'r')    
    lines = fl.readlines()
    
    # Remove comments
    for i, s in enumerate(lines):
        lines[i] = s.split("%", 1)[0]         
    
    # Remove lines with only whitespaces        
    lines = list(filter(lambda s: s.strip(), lines))        
    
    # Find all END keywords
    endIndices = [i for i, s in enumerate(lines) if 'end' in s.lower()]
    
    #===== Make a list of nodes (a Pandas datafram if things go right)
    nodes = []
    nodeBeginIndex = [i for i, s in enumerate(lines) if 'nodes' in s.lower()] # Find index where nodes list begins
    
    # Depending on the analysis, there won't be any node available
    if len(nodeBeginIndex) == 0:
        print(f"WARNING: Keyword 'nodes' missing from file {flPath}. Related analyzes, such as visualizing floater geommetry, will not be available.")        
        
    # Only one list of nodes can be included in the file
    elif len(nodeBeginIndex) > 1: 
        raise Exception(f"Keyword 'nodes' appearing more than once (not counting comments) in file {flPath}.")        
    
    else:
        # Jump index with only the 'NODES' keyword
        nodeBeginIndex = nodeBeginIndex[0]+1
        
        # There can be several END keywords in the file. We want the first one after NODES keyword
        nodeEndIndice = [i for i in endIndices if i > nodeBeginIndex]
        if len(nodeEndIndice) == 0:
            raise Exception(f"Need at least one end keyword after 'nodes' keyword in file {flPath}")
        nodeEndIndice = nodeEndIndice[0]
        
        # List of nodes stored as a dataframe
        nodes = pd.DataFrame([n.split(",") for n in lines[nodeBeginIndex:nodeEndIndice]])
        nodes.columns = ['id', 'x', 'y', 'z']
        
    #===== Make a list of circular Morison Elements (a Pandas dataframe if things go right)
    #===== The procedure is analogous to the one used for nodes
    morc = []
    morcBeginIndex = [i for i, s in enumerate(lines) if 'morison_circ' in s.lower()]
    
    if len(morcBeginIndex) == 0:
        print(f"WARNING: Keyword 'Morison_Circ' missing from file {flPath}. Related analyzes, such as visualizing floater geommetry, will not be available.")
        
    # Only one list of nodes can be included in the file
    elif len(morcBeginIndex) > 1:
        raise Exception(f"Keyword 'Morison_Circ' appearing more than once (not counting comments) in file {flPath}.")
    
    else:
        morcBeginIndex = morcBeginIndex[0]+1        
        morcEndIndice = [i for i in endIndices if i > morcBeginIndex]
        if len(morcEndIndice) == 0:
            raise Exception(f"Need at least one end keyword after 'Morison_Circ' keyword in file {flPath}")
        morcEndIndice = morcEndIndice[0]
        
        morc = pd.DataFrame([n.split() for n in lines[morcBeginIndex:morcEndIndice]])
        morc.columns = ['id1', 'id2', 'diam', 'cd', 'cm', 'npts', 'cdaxial1', 'caaxial1', 'cdaxial2', 'caaxial2', 'flagfk']
    
    
    
    # Dictionary with outputs
    out = {}
    out['nodes'] = nodes
    out['morc'] = morc
    return out
        

out = readInputFl('../test/OC4-JONSWAP.txt')

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
    


