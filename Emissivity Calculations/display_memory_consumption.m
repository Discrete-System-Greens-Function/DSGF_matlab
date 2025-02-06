function display_memory_consumption(vars)
% function to print out the workspace memory

	workspaceMem = sum(cell2mat(vars(3,:)))/1e9;
	disp('----------------------------------------------------------------')
	disp(['Workspace memory required is ', num2str(workspaceMem), ' GB'])
	disp('')
	disp('----------------------------------------------------------------')

end
