function [fowt, envir] = readInputFile(flPath)

    %= 'readInputFile' typed without arguments is the same as 'help readInputFile'
%     if nargin == 0
%       help(mfilename)
%       return
%     end

    if nargin == 0
        flPath = 'C:\Users\lucas.henrique\Google Drive\Doutorado\1Testes_Jappaku\bicromatica\Dir0_H2p0_BIC01-HD0.txt';
    end

    % open the input file
    fid = fopen(flPath, 'rt');
    if fid<0
       error('Couldn''t open file %s', flPath) 
    end
   
    lnC = {};
    envir.nodes.ID = [];
    envir.nodes.coord = [];
    fowt.floater.morisonElements = [];
    while ~feof(fid) % until the end of the file is achieved
        strInput  = readLineInputFile(fid);
        lnC{end+1} = strInput; % append current line to lnC

        if strcmpi(getKeyword(strInput), 'Nodes')
            % Read next line, since current line is just the main keyword
            strInput = readLineInputFile(fid);
    
            while ~strcmpi(getKeyword(strInput), 'END') 
                if (feof(fid))
                    error('End of file reached before END keyword in NODES specification.');
                end
    
                % Nodes are specified by a vec with four components: ID, X coord, Y coord, and Z coord.
                % They are separated by commas in the input string.
                inData = strsplit(strInput, ',');
                if numel(inData) ~= 4
                    error('Wrong number of parameters.');
                end
   
                envir.nodes.ID(end+1) = str2double(inData{1});
                envir.nodes.coord(end+1,:) = [str2double(inData{2}), str2double(inData{3}), str2double(inData{4})];
                
                % Done with this line. Read next one.
                strInput = readLineInputFile(fid);
            end        

        elseif strcmpi(getKeyword(strInput), 'FloaterCoG')       
            % The coordinates of the CoG are separated by commas in the input string
            inData = strsplit(getData(strInput), ',');
            if size(inData) ~= 3
                error('Unable to read the CoG. Wrong number of parameters.');
            end                
            fowt.floater.cog = [str2double(inData{1}), str2double(inData{2}), str2double(inData{3})];
   
        elseif strcmpi(getKeyword(strInput), 'Morison_circ')
            % Read next line, since current line is just the main keyword
            strInput = readLineInputFile(fid);
            
            while ~strcmpi(getKeyword(strInput), 'END') 
                if (feof(fid))
                    error('End of file reached before END keyword in MORISON_CIRC specification.');
                end
    
                % The eleven properties of a circular cylinder Morison's Element are separated by white spaces in the input string.
                inData = strsplit(strInput);
                if strcmp(inData(end), '')
                    inData(end) = [];
                end
                if numel(inData) ~= 11
                    error('Wrong number of parameters.');
                end
                
                if isempty(envir.nodes.ID)
                    error('Nodes should be specified before Morison Elements.');
                end
   
                indNode1 = find(envir.nodes.ID == str2double(inData{1}));
                indNode2 = find(envir.nodes.ID == str2double(inData{2}));
                
                if numel(indNode1) > 1
                    error(['Node ' num2str(envir.nodes.ID(indNode1(1))) ' was specified twice.' ])
                elseif numel(indNode2) > 1
                    error(['Node ' num2str(envir.nodes.ID(indNode2(1))) ' was specified twice.' ])
                end
                fowt.floater.morisonElements(end+1).node1 = envir.nodes.coord(indNode1,:);
                fowt.floater.morisonElements(end).node2 = envir.nodes.coord(indNode2,:);
                fowt.floater.morisonElements(end).diam = str2double(inData{3});
                fowt.floater.morisonElements(end).CD = str2double(inData{4});
                fowt.floater.morisonElements(end).CM = str2double(inData{5});
                fowt.floater.morisonElements(end).numIntPoints = str2double(inData{6});
                fowt.floater.morisonElements(end).botDiam = str2double(inData{7});
                fowt.floater.morisonElements(end).topDiam = str2double(inData{8});
                fowt.floater.morisonElements(end).axialCD = str2double(inData{9});
                fowt.floater.morisonElements(end).axialCa = str2double(inData{10});
                fowt.floater.morisonElements(end).botPressFlag = str2double(inData{11});
                                                 
                % Done with this line. Read next one.
                strInput = readLineInputFile(fid);                                                 
            end
            
%         elseif ~strcmpi(getKeyword(strInput), 'END_OF_INPUT_FILE')
%             warning(['Unknown keyword ' getKeyword(strInput) '.']);
        end
    end
end






%====
%==== Funcoes auxiliares
%====
function outStr = readLineInputFile(fid)
	outStr = '';

	% Read next file line to string strInput and update line number counter.
    % Repeat this process until the line has some content or end of file is achieved.
    while ~hasContent(outStr) && ~feof(fid)
		outStr = fgetl(fid);

		% Remove comments from line
        if thereIsCommentInString(outStr)
            outStr = removeComments(outStr);
        end
    end

    if feof(fid)	
		outStr = 'END_OF_INPUT_FILE';
    end
end

function keyword = getKeyword(inStr)
	% Get the first part of the string, the one before the first '\t' or white-space
    keyword = strtok(inStr);
end

function outStr = getData(inStr)
	% Check if input string is empty
    if isempty(inStr)
		error('Empty string passed to getData()');
    end    

	str_tokenized = strsplit(inStr);

	% If str_tokenized has only one element, i.e. only the keyword, then return an empty string
    if numel(str_tokenized) == 1
        outStr = str_tokenized{1};
		return
    end

	outStr = '';
    for ii = 2 : numel(str_tokenized)		
		outStr = [outStr , str_tokenized{ii}, sprintf('\t')];
    end
end

% Verify whether a string contains a comment, marked by a '%'
function bool = thereIsCommentInString(inStr)
    bool = contains(inStr, '%');
end

% Verify whether a string has content, i.e. if:
% 1) it is not empty
% 2) it is not just white spaces or tabs
% 3) it does not start with a comment mark ('%')
function bool = hasContent(inStr)	
    inStr = inStr(~isspace(inStr));

    if (isempty(inStr))
		bool = 0;
        return
    end
    
    if inStr(1) == '%'
        bool = 0;
        return
    end
    
    bool = 1;
end

function outStr = removeComments(inStr)    
	outStr = inStr(1:strfind(inStr, '%')-1);
end