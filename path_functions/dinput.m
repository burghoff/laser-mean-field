function h = dinput(fh,type,textlabel,stringval,callback)
% type: 'pushbutton' or 'edit'
% textlabel: for 'edit', a label to display above the input
% stringval: intial string value
% callback: callback function

figure(fh);
cntrlargs = {'Style',type,'String',stringval,'Callback',callback,'Tag','dinput'};

fkids = get(fh,'Children');
positions = [];
for ii=1:length(fkids)
    if isequal(fkids(ii).Tag,'dinput')
        positions(end+1,:) = fkids(ii).Position;
    end
end

mypos = [20,20,60,20];
if size(positions,1)>0
    [mv,ml]=max(positions(:,2));
    mypos(2) = mv + positions(ml,4)+3;
end

cntrlargs{end+1}='Position';
cntrlargs{end+1}=mypos;
h=uicontrol(cntrlargs{:});

if isequal(type,'edit')
    myposlabel = mypos;
    myposlabel(4)=15;
    myposlabel(2)=myposlabel(2)+20;
    uicontrol('Style','text','Position',myposlabel,'String',textlabel,'BackgroundColor',[1,1,1],'Tag','dinput');
end
end

