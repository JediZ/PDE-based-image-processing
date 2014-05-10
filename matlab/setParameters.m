%----------------------
%--- setParameters
%----------------------
function new_struct = setParameters(old_struct,varargin)
%function new_param = setParameters(old_param,varargin)
%
%Changes parameters in a parameter struct. The new parameters are given as
%['name',value] corresponding pairs in the varargin.
%
%INPUT
%old_struct     =       Old structure
%varargin       =       ['name',value]
%OUTPUT
%new_struct     =       New structure
%

n_param = length(varargin);
new_struct = old_struct;

for i=1:2:n_param

	%Extract parameter name
	p_name = varargin{i};

	%If for the parameter (name) exists a value...
	if i+1<=n_param
		%then extract the value
		p_value = varargin{i+1};
		
		%and if the p_name exists in the structure
		if isfield(new_struct, p_name)
			%set its value
			new_struct=setfield(new_struct, p_name, p_value);
		else
			%error(['setParameters error: unkown parameter name "',p_name,'". ',num2str(i)]);
		end
	end

end
