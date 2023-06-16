%fileName = 'F:\data\20210410 3D lipid+ATTO647N\';
fileName = 'E:\Experimental_data\20230208 SLB\';
ID = 33;
offsetName = [fileName,'_',num2str(ID),'\_',num2str(ID),'_MMStack_Default.ome.tif'];
Nimg = 100;

offsetR = Tiff(offsetName,'r');
for i=1:Nimg
    setDirectory(offsetR,i);
    offset(:,:,i) = double(offsetR.read);

end
offset_mean = mean(offset,3);
offset = offset_mean;
save([fileName 'processed data RoSEO\offSet.mat'],'offset')
%offset_ROI=offset_mean;imwrite(uint16(offset_ROI),'offset_ROI572_763.tif');


