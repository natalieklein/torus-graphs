function custom_node_placement(graph_handle,sub_reg_inds)
% This function adjusts the position of the 24 nodes to predefined locations
% input: 
% graph_handle is the plot handle of graph object
% sub_reg_inds has the indices corresponding to each subregion

%% customizing the location of the nodes in the network plot
graph_handle.XData(sub_reg_inds.ca3)=21:25;
graph_handle.YData(sub_reg_inds.ca3)=[5 4.7 5 4.7 5];
text(23,5.3,'CA3');axis off

graph_handle.XData(sub_reg_inds.dg)=2:9;
graph_handle.YData(sub_reg_inds.dg)=[6 5.7 6 5.7 6 5.7 6 5.7];
text(5,6.3,'DG');axis off

graph_handle.XData(sub_reg_inds.sub)=-2:0;
graph_handle.YData(sub_reg_inds.sub)=[4 3.7 4];
text(-7.5,3.8,'Sub');axis off

graph_handle.XData(sub_reg_inds.pfc)=12:19;
graph_handle.YData(sub_reg_inds.pfc)=[3 2.7 3 2.7 3 2.7 3 2.7];
text(15,2.5,'PFC');axis off
