function colorLabel(string)

%colorLabel - quick helper function that adds a label to the colorbar 

set(get(colorbar,'YLabel'),'String',string);
end