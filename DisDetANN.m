classdef DisDetANN < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        FileMenu                      matlab.ui.container.Menu
        HelpMenu                      matlab.ui.container.Menu
        AboutMenu                     matlab.ui.container.Menu
        InstructionMenu               matlab.ui.container.Menu
        TabGroup                      matlab.ui.container.TabGroup
        ParameterTab                  matlab.ui.container.Tab
        LogPanel                      matlab.ui.container.Panel
        TextArea                      matlab.ui.control.TextArea
        PlotframePanel                matlab.ui.container.Panel
        PlotButtonGroup               matlab.ui.container.ButtonGroup
        PlotXYZcoordinatesButton      matlab.ui.control.RadioButton
        PlotPointNormalButton         matlab.ui.control.RadioButton
        PlotPointCurvatureButton      matlab.ui.control.RadioButton
        PlotPointDensityButton        matlab.ui.control.RadioButton
        UIAxes                        matlab.ui.control.UIAxes
        ParameterPanel                matlab.ui.container.Panel
        PointDensityButton            matlab.ui.control.Button
        XYZcoordinatesButton          matlab.ui.control.Button
        PointnormalButton             matlab.ui.control.StateButton
        PointcurvatureButton          matlab.ui.control.StateButton
        kEditFieldLabel               matlab.ui.control.Label
        kEditField                    matlab.ui.control.NumericEditField
        radiusEditFieldLabel          matlab.ui.control.Label
        radiusEditField               matlab.ui.control.NumericEditField
        MeandistanceButton            matlab.ui.control.Button
        meandistacneEditField         matlab.ui.control.NumericEditField
        LoadingPanel                  matlab.ui.container.Panel
        ViewrawdataButton             matlab.ui.control.Button
        FilepathTextAreaLabel         matlab.ui.control.Label
        FilepathTextArea              matlab.ui.control.TextArea
        NetworkTab                    matlab.ui.container.Tab
        logPanel                      matlab.ui.container.Panel
        TextArea_2                    matlab.ui.control.TextArea
        TraintheNetworkPanel          matlab.ui.container.Panel
        Panel_2                       matlab.ui.container.Panel
        TrainButton                   matlab.ui.control.StateButton
        SavenetButton                 matlab.ui.control.Button
        TrainFunDropDownLabel         matlab.ui.control.Label
        TrainFunDropDown              matlab.ui.control.DropDown
        PerformFunDropDownLabel       matlab.ui.control.Label
        PerformFunDropDown            matlab.ui.control.DropDown
        HiddenSizesEditFieldLabel     matlab.ui.control.Label
        HiddenSizesEditField          matlab.ui.control.NumericEditField
        Panel_3                       matlab.ui.container.Panel
        LoadnetButton                 matlab.ui.control.Button
        Panel_4                       matlab.ui.container.Panel
        PredictionButton              matlab.ui.control.StateButton
        SelecttheSamplesPanel         matlab.ui.container.Panel
        GroupsinthisoutcropEditFieldLabel  matlab.ui.control.Label
        GroupsinthisoutcropEditField  matlab.ui.control.NumericEditField
        GroupDropDownLabel            matlab.ui.control.Label
        GroupDropDown                 matlab.ui.control.DropDown
        SelectButton                  matlab.ui.control.Button
        PlotframePanel_2              matlab.ui.container.Panel
        SamplesButton_2               matlab.ui.control.Button
        ResultButton_4                matlab.ui.control.Button
        UIAxes2                       matlab.ui.control.UIAxes
        NormalizationPanel            matlab.ui.container.Panel
        normaliztionButton            matlab.ui.control.Button
        ClusteringTab                 matlab.ui.container.Tab
        ClusteringPanel               matlab.ui.container.Panel
        GroupDropDown_2Label          matlab.ui.control.Label
        GroupDropDown_2               matlab.ui.control.DropDown
        minPtsEditFieldLabel          matlab.ui.control.Label
        minPtsEditField               matlab.ui.control.NumericEditField
        Label                         matlab.ui.control.Label
        epsEditField                  matlab.ui.control.NumericEditField
        ClusterButton                 matlab.ui.control.Button
        ResultButton_2                matlab.ui.control.Button
        LogPanel_2                    matlab.ui.container.Panel
        TextArea_3                    matlab.ui.control.TextArea
        PotframePanel                 matlab.ui.container.Panel
        UIAxes3                       matlab.ui.control.UIAxes
        OrientationTab                matlab.ui.container.Tab
        ExportButton                  matlab.ui.control.Button
        ResultButton_3                matlab.ui.control.Button
        CalculationButton             matlab.ui.control.Button
        StereographicProjectionPanel  matlab.ui.container.Panel
        UIAxes4                       matlab.ui.control.UIAxes
        LogPanel_3                    matlab.ui.container.Panel
        TextArea_4                    matlab.ui.control.TextArea
        GroupDropDown_3               matlab.ui.control.DropDown
        GroupDropDown_3Label          matlab.ui.control.Label
    end

    
    properties (Access = private)
        pcData % Description
        % Description
        StrArray % Description
        ptCloud % Description
        pcDataNew % Description
        pcLearn % cell.
        net % Description
        outputs_test % Description
        groupNum % Description
        jointsData % Description
        Orientation % Description
    end
    
    methods (Access = private)
        
        function logRefresh_func(app,StrArrayNew)
            if length(app.StrArray)>=20
                app.StrArray={};
            end
            selectedTab=app.TabGroup.SelectedTab;
            
            switch selectedTab
                case app.ParameterTab
                    app.StrArray=[app.StrArray,StrArrayNew];
                    app.TextArea.Value=app.StrArray;
                case app.NetworkTab
                    app.StrArray=[app.StrArray,StrArrayNew];
                    app.TextArea_2.Value=app.StrArray;
                case app.ClusteringTab
                    app.StrArray=[app.StrArray,StrArrayNew];
                    app.TextArea_3.Value=app.StrArray;
                case app.OrientationTab
                    app.StrArray=[app.StrArray,StrArrayNew];
                    app.TextArea_4.Value=app.StrArray;
            end
            return
        end
        
        
        function  [s,v]=CovarianceMatrix(app,ptCloud,n)
            m=ptCloud.Count;
            h=waitbar(0,['Point curvature calculation in process. ',num2str(m),' points. Please wait']);
            v=zeros(3,ptCloud.Count);
            s=zeros(ptCloud.Count,1);
            % p=zeros(3,ptCloud.Count);
            for i=1:ptCloud.Count
                [indices,~] = findNearestNeighbors(ptCloud,ptCloud.Location(i,:),n);
                x = ptCloud.Location(indices(:),:);
                p_bar = 1/n * sum(x,1);
                P = transpose(x - repmat(p_bar,n,1))*(x - repmat(p_bar,n,1)); 
                [V,lmd] = eig(P);
                [lmds,id]=min(diag(lmd));
                v(:,i)=V(:,id);
                s(i) = lmds/sum(diag(lmd))*100;
                waitbar(i/m,h);
            end
            
            close(h);
        end
        
        
        function pos = getpointsXYZ(app,data,n)
            h = figure;
            pcshow(data(:,1:3),data(:,4:6));
            grid on;
            set(gca,'fontname','Times New Roman','fontsize',14);
            xlabel(gca,'X (m)','fontname','Times New Roman','fontsize',16 );
            ylabel(gca,'Y (m)','fontname','Times New Roman','fontsize',16 );
            zlabel(gca,'Z (m)','fontname','Times New Roman','fontsize',16 );
            axis equal;
            view(340,10);
            datacursormode on
            dcm_obj = datacursormode(h);
            pos = zeros(n,3);
            for i = 1:n
                w = waitforbuttonpress;
                while w==0
                    w = waitforbuttonpress;
                end
                c_info = getCursorInfo(dcm_obj);
                pos(i,:) = c_info.Position;
            end
            close(h);
        end
        
         
        
        function [ T ] = f_dbscan( app,A , eps, ppcluster)%
            % [ T, eps ] = f_dbscan( A , npb, ppcluster)
            % Búsqueda de clústers mediante una búsqueda previa de vecinos
            % Aplicación del algoritmo DBSCAN
            % Adrián Riquelme Guill, mayo 2013
            %    Copyright (C) {2015}  {Adrián Riquelme Guill, adririquelme@gmail.com
            
            [n,d]=size(A);
            h=waitbar(0,['Cluster analysis in process. ',num2str(n),' points. Please wait']);
            
            minpts=d+1; %  
            T=zeros(n,1);
            maxcluster=1; % 
            [idx, ~] = rangesearch(A,A,eps);
            for i=1:n
                NeighborPts=idx{i};
                if length(NeighborPts)>=minpts 
                    cv=T(NeighborPts);
                    mincv=min(cv); 
                    mincv2=min(cv((cv>0))); 
                    maxcv=max(cv);
                    if maxcv==0
                        caso=0;  
                    else
                        if maxcv==mincv2
                            caso=1; 
                        else
                            caso=2; 
                        end
                    end
                    switch caso
                        case 0
                            % 
                            T(NeighborPts)=maxcluster; 
                            % T(i)=maxcluster;
                            maxcluster=maxcluster+1; %
                        case 1
                            if mincv==0
                                T(NeighborPts(cv==0))=mincv2;
                            end
                            % T(i)=mincv2;
                        case 2
                            T(NeighborPts(cv==0))=mincv2;
                            % 
                            b=cv(cv>mincv2); 
                            [~,n1]=size(b);
                            aux=0;
                            for j=1:n1
                                if b(j)~=aux
                                    T(T==b(j))=mincv2;
                                    aux=b(j);
                                end
                            end
                            % T(i)=mincv2;
                    end
                else
                end
                waitbar(i/n,h);
            end
            %% 

            if sum(T)==0

            else
               
                T2=T;
                cluster=unique(T2,'sorted');
                cluster=cluster(cluster>0); 
                [ nclusters,~]=size(cluster);
                A=zeros(2,nclusters);
                numeroclusters=zeros(1, nclusters);
                for ii=1:nclusters
                    numeroclusters(ii)=length(find(T2(:,1)==cluster(ii,1)));
                end
                A(2,:)=cluster; A(1,:)=numeroclusters;   
               
                [~,IX]=sort(A(1,:),'descend'); A=A(:,IX);

                n=ppcluster;
                I=find(A(1,:)>n);
                J=find(A(1,:)<=n);

                for ii=1:length(J)
                    T(T2==A(2,J(ii)))=0;
                end

                for ii=1:length(I)
                    T(T2==A(2,I(ii)))=ii;
                end
            end
            close(h);
        end
        
        function Vector=PointCloudVector(app,points)

 
            b=zeros(1,3)-sum(points,1);
            a=points'*points;
            N=a\b';

            Vector=(N/norm(N))';
        end
        
        function [dip,dd] = OrientationM(app,i)

            x=i(:,1);
            y=i(:,2);
            z=i(:,3);
            xy=sqrt(x.^2+y.^2);
            
            if z<0
                dip=pi-acos(z);
                if x>=0 && y>=0
                    dd=asin(x/xy)+pi;
                elseif x<0 && y>=0
                    dd=2*pi-asin(-x/xy)-pi;
                elseif x>=0 && y<0
                    dd=pi-asin(x/xy)+pi;
                else
                    %elseif x<0 && y<0
                    dd=pi+asin(-x/xy)-pi;
                end
            else
                dip=acos(z);
                if x>=0 && y>=0
                    dd=asin(x/xy);
                elseif x<0 && y>=0
                    dd=2*pi-asin(-x/xy);
                elseif x>=0 && y<0
                    dd=pi-asin(x/xy);
                else
                    %elseif x<0 && y<0
                    dd=pi+asin(-x/xy);
                end
            end
            
            dip=rad2deg(dip);
            dd=rad2deg(dd);
            
        end
        
        function Stereonet( app,ori,r,int,colorindex )
            
            i=0:int:360;
            j=0:int:90;
            [i,j]=meshgrid(i,j);
            x=r.*sin(deg2rad(i)).*tan(deg2rad(j)/2);
            y=r.*cos(deg2rad(i)).*tan(deg2rad(j)/2);
            z = zeros(size(x));
            mesh(app.UIAxes4,y,z,'EdgeColor','k','FaceAlpha',0,'LineWidth',0.5);
            hold (app.UIAxes4,"on");
            grid (app.UIAxes4,"off");
            view(app.UIAxes4,2)
            axis (app.UIAxes4,"equal");
            
            
            dip=ori(:,1);
            dd=ori(:,2);
            xx=r.*sin(deg2rad(dd)).*tan(deg2rad(dip)/2);
            yy=r.*cos(deg2rad(dd)).*tan(deg2rad(dip)/2);
            plot(app.UIAxes4,xx,yy,'o','MarkerSize',8,'MarkerFaceColor',colorindex,'MarkerEdgeColor',colorindex);
            text(app.UIAxes4,-5,110,'N','FontSize',14);
            text(app.UIAxes4,-5,-110,'S','FontSize',14);
            text(app.UIAxes4,-120,0,'W','FontSize',14);
            text(app.UIAxes4,105,0,'E','FontSize',14);
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.UIFigure.Name='Discontinuity Detection using ANN (DisDetANN) V 1.0';
        end

        % Button pushed function: ViewrawdataButton
        function ViewrawdataButtonPushed(app, event)
            scatter3(app.UIAxes,app.pcData(:,1),app.pcData(:,2),app.pcData(:,3),3,app.pcData(:,3),'filled');
            grid (app.UIAxes,"on");
            set(app.UIAxes,'fontname','Times New Roman','fontsize',14);
            xlabel(app.UIAxes,'X (m)','fontname','Times New Roman','fontsize',16 );
            ylabel(app.UIAxes,'Y (m)','fontname','Times New Roman','fontsize',16 );
            zlabel(app.UIAxes,'Z (m)','fontname','Times New Roman','fontsize',16 );
            axis (app.UIAxes,"equal");
            view(app.UIAxes,45,10);
            title(app.UIAxes,'Raw data')
            title(app.UIAxes,'Raw data')
        end

        % Button pushed function: XYZcoordinatesButton
        function XYZcoordinatesButtonPushed(app, event)
            app.pcData=app.pcData;
            StrArrayNew={'Calculation is completed！'};
            logRefresh_func(app,StrArrayNew)           
        end

        % Value changed function: PointnormalButton
        function PointnormalButtonValueChanged(app, event)
            app.ptCloud=pointCloud(app.pcData(:,1:3));
            k=app.kEditField.Value;
            pitNormal=pcnormals(app.ptCloud,k);
            app.pcData(:,4:6)=pitNormal;
            StrArrayNew={'Point normal calculation is completed！'};
            logRefresh_func(app,StrArrayNew)
        end

        % Value changed function: PointcurvatureButton
        function PointcurvatureButtonValueChanged(app, event)
            k=app.kEditField.Value;
            [s,~]=CovarianceMatrix(app,app.ptCloud,k);
            app.pcData(:,7)=s;
            StrArrayNew={'Point curvature calculation is completed！'};
            logRefresh_func(app,StrArrayNew)
        end

        % Button pushed function: PointDensityButton
        function PointDensityButtonPushed(app, event)
            K=app.radiusEditField.Value;
            n=size(app.ptCloud.Location,1);
            h=waitbar(0,['Point density calculation in process. ',num2str(n),' points. Please wait']);
            pDen=zeros(n,1);
            for i=1:n
                [indices,~] = findNeighborsInRadius(app.ptCloud,app.ptCloud.Location(i,:),K);
                pDen(i,1)=size(indices,1);
                waitbar(i/n,h);
            end
            close(h);
            app.pcData(:,8)=pDen;
            StrArrayNew={'Point density calculation is completed！'};
            logRefresh_func(app,StrArrayNew)
        end

        % Selection changed function: PlotButtonGroup
        function PlotButtonGroupSelectionChanged(app, event)

            selectedButton = app.PlotButtonGroup.SelectedObject;
            switch selectedButton.Text
                case 'XYZ-coordinates'
                    scatter3(app.UIAxes,app.pcData(:,1),app.pcData(:,2),app.pcData(:,3),3,app.pcData(:,1:3),'filled');
                    grid (app.UIAxes,"on");
                    set(app.UIAxes,'fontname','Times New Roman','fontsize',14);
                    xlabel(app.UIAxes,'X (m)','fontname','Times New Roman','fontsize',16 );
                    ylabel(app.UIAxes,'Y (m)','fontname','Times New Roman','fontsize',16 );
                    zlabel(app.UIAxes,'Z (m)','fontname','Times New Roman','fontsize',16 );
                    axis (app.UIAxes,"equal");
                    view(app.UIAxes,45,10);
                    title(app.UIAxes,'Raw data with XYZ-coordinate as point cloud color')
                    StrArrayNew={'Image displayed and the color is the original coordinate！'};
                    logRefresh_func(app,StrArrayNew)
                case 'Point Normal'
                    scatter3(app.UIAxes,app.pcData(:,1),app.pcData(:,2),app.pcData(:,3),3,app.pcData(:,4:6),'filled');
                    grid (app.UIAxes,"on");
                    set(app.UIAxes,'fontname','Times New Roman','fontsize',14);
                    xlabel(app.UIAxes,'X (m)','fontname','Times New Roman','fontsize',16 );
                    ylabel(app.UIAxes,'Y (m)','fontname','Times New Roman','fontsize',16 );
                    zlabel(app.UIAxes,'Z (m)','fontname','Times New Roman','fontsize',16 );
                    axis (app.UIAxes,"equal");
                    view(app.UIAxes,45,10);
                    title(app.UIAxes,'Raw data with Point Normal as point cloud color')
                    StrArrayNew={'Image displayed and the color is the point normal！'};
                    logRefresh_func(app,StrArrayNew)
                case 'Point Curvature'
                    scatter3(app.UIAxes,app.pcData(:,1),app.pcData(:,2),app.pcData(:,3),3,app.pcData(:,7),'filled');
                    grid (app.UIAxes,"on");
                    set(app.UIAxes,'fontname','Times New Roman','fontsize',14);
                    xlabel(app.UIAxes,'X (m)','fontname','Times New Roman','fontsize',16 );
                    ylabel(app.UIAxes,'Y (m)','fontname','Times New Roman','fontsize',16 );
                    zlabel(app.UIAxes,'Z (m)','fontname','Times New Roman','fontsize',16 );
                    axis (app.UIAxes,"equal");
                    view(app.UIAxes,45,10);
                    title(app.UIAxes,'Raw data with Point Curvature as point cloud color')
                    StrArrayNew={'Image displayed and the color is the point curvature！'};
                    logRefresh_func(app,StrArrayNew)
                case 'Point Density'
                    scatter3(app.UIAxes,app.pcData(:,1),app.pcData(:,2),app.pcData(:,3),3,app.pcData(:,8),'filled');
                    grid (app.UIAxes,"on");
                    set(app.UIAxes,'fontname','Times New Roman','fontsize',14);
                    xlabel(app.UIAxes,'X (m)','fontname','Times New Roman','fontsize',16 );
                    ylabel(app.UIAxes,'Y (m)','fontname','Times New Roman','fontsize',16 );
                    zlabel(app.UIAxes,'Z (m)','fontname','Times New Roman','fontsize',16 );
                    axis (app.UIAxes,"equal");
                    view(app.UIAxes,45,10);
                    title(app.UIAxes,'Raw data with Point Density as point cloud color')
                    StrArrayNew={'Image displayed and the color is the point density！'};
                    logRefresh_func(app,StrArrayNew)
            end
            
        end

        % Button pushed function: normaliztionButton
        function normaliztionButtonPushed(app, event)
            [m,n]=size(app.pcData);
            app.pcDataNew=zeros(m,n);
            for i=1:1:n
                maxVar=max(app.pcData(:,i));
                minVar=min(app.pcData(:,i));
                for j=1:1:m
                    app.pcDataNew(j,i)=(app.pcData(j,i)-minVar)/(maxVar-minVar);
                end
            end
            StrArrayNew={'Data normalization is completed！'};
            logRefresh_func(app,StrArrayNew)
        end

        % Value changed function: GroupsinthisoutcropEditField
        function GroupsinthisoutcropEditFieldValueChanged(app, event)
            app.groupNum=app.GroupsinthisoutcropEditField.Value;
            app.Orientation = cell(app.groupNum-1,1);
            for i=1:app.groupNum
                if i==1
                    NameofGroup='Non-joint';
                    app.GroupDropDown.Items{1,i}=NameofGroup;
                else
                    th=num2str(i-1);
                    NameofGroup=strcat('Joint set',32,th);
                    app.GroupDropDown.Items{1,i}=NameofGroup;
                    app.GroupDropDown_2.Items{1,i-1}=NameofGroup;
                    app.GroupDropDown_3.Items{1,i-1}=NameofGroup;    
                end
            end
            app.pcLearn=cell(1,app.groupNum);
        end

        % Button pushed function: SelectButton
        function SelectButtonPushed(app, event)
            value = app.GroupDropDown.Value;
            num_str = regexp(value,'\d*\.?\d*','match');
            num = str2double(num_str);
            prompt = {'How many points do you want to select for this data group? '};
            dlgtitle = 'Input';
            dims = [1 50];
            answer = inputdlg(prompt,dlgtitle,dims);
            pointsNum= str2num(answer{1}); %#ok<*ST2NM> 
            pointsGet=getpointsXYZ(app,app.pcData,pointsNum);           
            pointsGet(:,10:9+app.groupNum)=zeros(size(pointsGet,1),app.groupNum);
            if isempty(num)
                num=0;
                pointsGet(:,10+num)=ones(pointsNum,1);%
            else
                pointsGet(:,10+num)=ones(pointsNum,1);%
            end
            [~,Locb]=ismember(pointsGet(:,1:3),app.pcData(:,1:3),'rows');
            pointsGet(:,1:8)=app.pcDataNew(Locb,1:8);
            pointsGet(:,9)=Locb;
            app.pcLearn{num+1}=pointsGet;
            StrArrayNew=strcat(num2str(pointsNum),' points were selected for',32,value);
            logRefresh_func(app,StrArrayNew);
            StrArrayNew=strcat('Selecting data in the',32, value,32 ,'is completed');
            logRefresh_func(app,StrArrayNew)
        end

        % Value changed function: TrainButton
        function TrainButtonValueChanged(app, event)
            pclearn=cat(1,app.pcLearn{:});
            inputs = pclearn(:,1:8)';
            targets = pclearn(:,10:9+app.groupNum)';
            
            % Create a Pattern Recognition Network
            hiddenSizes =app.HiddenSizesEditField.Value;
            trainFcn=app.TrainFunDropDown.Value;
            performFcn=app.PerformFunDropDown.Value;            
            app.net=patternnet(hiddenSizes,trainFcn,performFcn);         
            % Train the Network
            [app.net,~] = train(app.net,inputs,targets);
            nntraintool;
            nntraintool('close');
            StrArrayNew=('Network training has been completed');
            logRefresh_func(app,StrArrayNew)
        end

        % Value changed function: PredictionButton
        function PredictionButtonValueChanged(app, event)
          app.outputs_test = round(app.net(app.pcDataNew(:,1:8)'))';
          app.pcData(:,1:8+app.groupNum)= [app.pcData(:,1:8) app.outputs_test];
          StrArrayNew={'Classfication is completed'};
          logRefresh_func(app,StrArrayNew)
        end

        % Button pushed function: ClusterButton
        function ClusterButtonPushed(app, event)
            eps=app.epsEditField.Value;
            ppcluster=app.minPtsEditField.Value;
            value = app.GroupDropDown_2.Value;
            num_str = regexp(value,'\d*\.?\d*','match');
            num = str2double(num_str);
            if isempty(app.outputs_test)
            else
            jointset=find(app.outputs_test(:,num+1)==1);
            end
            if isempty(app.pcData)
            else            
            app.pcData(jointset,9+app.groupNum)=f_dbscan(app,app.pcData(jointset,1:3),eps,ppcluster);
            end
            StrArrayNew=strcat('The clustering of the',32,value,32,'has been completed');
            logRefresh_func(app,StrArrayNew)
        end

        % Button pushed function: ResultButton_2
        function ResultButton_2Pushed(app, event)
            for ii=1:1:app.groupNum-1
                jointset=find(app.pcData(:,9+ii)==1);
                m=max(app.pcData(jointset,9+app.groupNum));
                cx=rand(m,1);
                cy=rand(m,1);
                cz=rand(m,1);
                for jj=1:m
                    color(:,1)= cx(jj);
                    color(:,2)= cy(jj);
                    color(:,3)= cz(jj);
                    j=app.pcData(find(app.pcData(:,9+ii)==1&app.pcData(:,9+app.groupNum)==jj),:);
                    scatter3(app.UIAxes3,j(:,1),j(:,2),j(:,3),3,color,'filled');
                    hold (app.UIAxes3,"on");
                end 
            end
            grid (app.UIAxes3,"on");
            set(app.UIAxes3,'fontname','Times New Roman','fontsize',14);
            xlabel(app.UIAxes3,'X (m)','fontname','Times New Roman','fontsize',16 );
            ylabel(app.UIAxes3,'Y (m)','fontname','Times New Roman','fontsize',16 );
            zlabel(app.UIAxes3,'Z (m)','fontname','Times New Roman','fontsize',16 );
            axis (app.UIAxes3,"equal");
            view(app.UIAxes3,45,10);
            title(app.UIAxes3,'Clustering results')
            hold (app.UIAxes3,"off");
            StrArrayNew=('Clustering results have been shown and different individual joint has different color ');
            logRefresh_func(app,StrArrayNew)
        end

        % Button pushed function: CalculationButton
        function CalculationButtonPushed(app, event)
            value = app.GroupDropDown_3.Value;
            num_str = regexp(value,'\d*\.?\d*','match');
            num = str2double(num_str);
            jointset=find(app.pcData(:,9+num)==1);
            m=max(app.pcData(jointset,9+app.groupNum));
            for jj=1:m
                j=app.pcData(find(app.pcData(:,9+num)==1&app.pcData(:,9+app.groupNum)==jj),:);
                jointNor=PointCloudVector(app,j(:,1:3));
                [dip,dipdd]=OrientationM(app,jointNor);
                [n,~]=size(j);
                app.Orientation{num}(jj,1)=jj;
                app.Orientation{num}(jj,2)=n;
                app.Orientation{num}(jj,3)=dip;
                app.Orientation{num}(jj,4)=dipdd;
            end
            StrArrayNew=strcat('The orientation calculation of',32,value,32,'is completed');
            logRefresh_func(app,StrArrayNew)
        end

        % Button pushed function: MeandistanceButton
        function MeandistanceButtonPushed(app, event)
            % Búsqueda de clústers mediante una búsqueda previa de vecinos
            % Aplicación del algoritmo DBSCAN
            % Adrián Riquelme Guill, mayo 2013
            %    Copyright (C) {2015}  {Adrián Riquelme Guill, adririquelme@gmail.com
            ksigmaseps = 2;
            nvecinos=app.kEditField.Value+1;
            [n,~]=size(app.pcData(:,1:3));
            if nvecinos > n  
                nvecinos=n;
                [~,dist]=knnsearch(app.pcData(:,1:3),app.pcData(:,1:3),'NSMethod','kdtree','distance','euclidean','k',nvecinos);
                data=dist(:,nvecinos); 
            else
                [~,dist]=knnsearch(app.pcData(:,1:3),app.pcData(:,1:3),'NSMethod','kdtree','distance','euclidean','k',nvecinos);
                if n<5
                    data=dist(:,n); 
                else                                                               
                    data=dist(:,5); %
                end
            end
            data=unique(data,'sorted'); 
            Meandistance=mean(data)+ksigmaseps*std(data);
            Meandistance=double(Meandistance);
            app.meandistacneEditField.Value=Meandistance;
            StrArrayNew={'Mean distance calculation is completed！'};
            logRefresh_func(app,StrArrayNew)
        end

        % Button pushed function: ResultButton_3
        function ResultButton_3Pushed(app, event)
            r=100;int=15;
            i=0:int:360;
            j=0:int:90;
            [i,j]=meshgrid(i,j);
            x=r.*sin(deg2rad(i)).*tan(deg2rad(j)/2);
            y=r.*cos(deg2rad(i)).*tan(deg2rad(j)/2);
            z = zeros(size(x));
            mesh(app.UIAxes4,x,y,z,'EdgeColor','k','FaceAlpha',0,'LineWidth',0.5);
            hold (app.UIAxes4,"on");
            grid (app.UIAxes4,"off");
            axis(app.UIAxes4,"off")
            view(app.UIAxes4,2)
            axis (app.UIAxes4,"equal");
            cx=rand(app.groupNum-1,1);
            cy=rand(app.groupNum-1,1);
            cz=rand(app.groupNum-1,1);
            for ii=1:app.groupNum-1
                ori=app.Orientation{ii}(:,3:4);
                color(:,1)= cx(ii);
                color(:,2)= cy(ii);
                color(:,3)= cz(ii);
                dip=ori(:,1);
                dd=ori(:,2);
                xx=r.*sin(deg2rad(dd)).*tan(deg2rad(dip)/2);
                yy=r.*cos(deg2rad(dd)).*tan(deg2rad(dip)/2);
                plot(app.UIAxes4,xx,yy,'o','MarkerSize',8,'MarkerFaceColor',color,'MarkerEdgeColor',color);
                hold(app.UIAxes4,"on")
            end
            text(app.UIAxes4,-5,110,'N','FontSize',14);
            text(app.UIAxes4,-5,-110,'S','FontSize',14);
            text(app.UIAxes4,-120,0,'W','FontSize',14);
            text(app.UIAxes4,105,0,'E','FontSize',14);
            title(app.UIAxes4,'');
            hold(app.UIAxes4,"off")
        end

        % Button pushed function: SavenetButton
        function SavenetButtonPushed(app, event)
            [file,path] = uiputfile({'*.mat','mat-files(*.mat)'},'save');
            str=fullfile(path,file);
            temp=app.net;
            save(str,'temp');             
        end

        % Button pushed function: LoadnetButton
        function LoadnetButtonPushed(app, event)
            [file,path] = uigetfile({'*.mat','mat-files(*.mat)'},'load');
            str=fullfile(path,file);
            S=load(str);
            names = fieldnames(S);
            app.net= S.(names{1});  
            StrArrayNew={'The network has been imported！'};
            logRefresh_func(app,StrArrayNew)
        end

        % Button pushed function: ExportButton
        function ExportButtonPushed(app, event)
            numJoint=0;
            for ii=1:app.groupNum-1
                [n,~]=size(app.Orientation{ii});
                for jj=1:n
                app.Orientation{ii}(jj,1)=jj;
                end
                numJoint=numJoint+n;
            end
            Orientation_temp=cat(1,app.Orientation{:});
            T=table(Orientation_temp(:,1),Orientation_temp(:,2),Orientation_temp(:,3),Orientation_temp(:,4),'VariableNames',{'Num','Count','Dip','Dip Direction'});
            [file,path] = uiputfile({'*.xlsx','xlsx-files(*.xlsx)'},'save');
            str=fullfile(path,file);
            writetable(T,str);
        end

        % Menu selected function: FileMenu
        function FileMenuSelected(app, event)
            [FileName,PathName,~] = uigetfile({'*.xls';'*.xlsx';'*.txt';'*.mat';'*.pcd'},'Please selectpoint cloud data');
            app.FilepathTextArea.Value=fullfile(PathName,FileName);
            switch(FileName(end-2:end))  
                case 'mat'
                    data=load([PathName,FileName]);
                    names = fieldnames(data);
                    app.pcData= data.(names{1}); 
                    StrArrayNew={'Data imported successfully！'};
                    logRefresh_func(app,StrArrayNew)
                case 'xls'
                    app.pcData=xlsread([PathName,FileName]);
                    StrArrayNew={'Data imported successfully！'};
                    logRefresh_func(app,StrArrayNew)
                case 'lsx'
                    app.pcData=xlsread([PathName,FileName]);
                    StrArrayNew={'Data imported successfully！'};
                    logRefresh_func(app,StrArrayNew)
                case 'txt'
                    data=importdata([PathName,FileName]);
                    app.pcData=data(:,1:3);
                    StrArrayNew={'Data imported successfully！'};
                    logRefresh_func(app,StrArrayNew)
                case 'pcd'
                    data_temp=pcread([PathName,FileName]);
                    app.pcData=data_temp.Location;
                    StrArrayNew={'Data imported successfully！'};
                    logRefresh_func(app,StrArrayNew)
                otherwise
                    StrArrayNew={'Only.mat/.xls/.xlsx/.pcd data is supported. This type is not currently supported！'};
                    logRefresh_func(app,StrArrayNew)
            end
            return
        end

        % Menu selected function: AboutMenu
        function AboutMenuSelected(app, event)


    d = dialog('Position',[300 300 800 400],'Name','About');

    uicontrol('Parent',d,...
               'Style','text', ...
               'FontName','Times New Roman',...
               'String', 'Discontinuity Detection using ANN with 3D point clouds, 2021', ...   
               'Position',[10 280 800 20],...
               'FontSize',12);

    uicontrol('Parent',d,...
               'Style','text',...
               'FontName','Times New Roman',...
               'Position',[10 255 800 20],...
               'FontSize',12,...
               'String','Yunfeng Ge, Bei Cao. China University of Geosciences');
    
    uicontrol('Parent',d,...
               'Style','text',...
               'FontName','Times New Roman',...
               'Position',[10 230 800 20],...
               'FontSize',12,...
               'String','geyunfeng@cug.edu.cn');
    
    uicontrol('Parent',d,...
               'Style','text',...
               'FontName','Times New Roman',...
               'Position',[10 190 800 20],...
               'FontSize',12,...
               'String','Cite as:');
    
    uicontrol('Parent',d,...
               'Style','text',...
               'FontName','Times New Roman',...
               'Position',[10 165 800 20],...
               'FontSize',12,...
               'String','Yunfeng Ge, Bei Cao, Huiming Tang, Binbin Zhao. Automated identification of rock discontinuities ');
    
    uicontrol('Parent',d,...
               'Style','text',...
               'FontName','Times New Roman',...
               'Position',[10 140 800 20],...
               'FontSize',12,...
               'String','from 3D point clouds using atificial intelligence. Manuscript Number: RMRE-D-21-00558');
    
     uicontrol('Parent',d,...
               'Style','text',...
               'FontName','Times New Roman',...
               'Position',[10 115 800 20],...
               'FontSize',12,...
               'String','Keywords: rock discontinuities; automated identification; point clouds; artificail nerual network');




        end

        % Button pushed function: SamplesButton_2
        function SamplesButton_2Pushed(app, event)
            cx=rand(app.groupNum,1);
            cy=rand(app.groupNum,1);
            cz=rand(app.groupNum,1);
            
            for i=1:app.groupNum
                pointselected_temp=app.pcLearn{i};
                if isempty(pointselected_temp)
                    
                else
                    [~,b]=ismember(pointselected_temp(:,1:3),app.pcDataNew(:,1:3),'rows');
                    pointselected=app.pcData(b,:);
                    n=size(pointselected,1);
                    color=zeros(n,3);
                    color(:,1)= cx(i);
                    color(:,2)= cy(i);
                    color(:,3)= cz(i);
                    scatter3(app.UIAxes2,pointselected(:,1),pointselected(:,2),pointselected(:,3),40,color,'filled');
                    hold (app.UIAxes2,"on");
                end
            end
           
            scatter3(app.UIAxes2,app.pcData(:,1),app.pcData(:,2),app.pcData(:,3),3,'filled');
            grid (app.UIAxes2,"on");
            set(app.UIAxes2,'fontname','Times New Roman','fontsize',14);
            xlabel(app.UIAxes2,'X (m)','fontname','Times New Roman','fontsize',16 );
            ylabel(app.UIAxes2,'Y (m)','fontname','Times New Roman','fontsize',16 );
            zlabel(app.UIAxes2,'Z (m)','fontname','Times New Roman','fontsize',16 );
            axis (app.UIAxes2,"equal");
            view(app.UIAxes2,45,10);
            title(app.UIAxes2,'Samples')
            hold (app.UIAxes2,"off");
            StrArrayNew={'Samples displayed and the points belonging to the same group were assigned to same color！'};
            logRefresh_func(app,StrArrayNew)
            
        end

        % Button pushed function: ResultButton_4
        function ResultButton_4Pushed(app, event)
            cx=rand(app.groupNum,1);
            cy=rand(app.groupNum,1);
            cz=rand(app.groupNum,1);
            outputs_color=zeros(size(app.outputs_test,1),3);
            for i=1:app.groupNum
                outputs_color(app.outputs_test(:,i)==1,1) = cx(i);
                outputs_color(app.outputs_test(:,i)==1,2) = cy(i);
                outputs_color(app.outputs_test(:,i)==1,3) = cz(i);
            end
            scatter3(app.UIAxes2,app.pcData(:,1),app.pcData(:,2),app.pcData(:,3),3,outputs_color,'filled');
            grid (app.UIAxes2,"on");
            set(app.UIAxes2,'fontname','Times New Roman','fontsize',14);
            xlabel(app.UIAxes2,'X (m)','fontname','Times New Roman','fontsize',16 );
            ylabel(app.UIAxes2,'Y (m)','fontname','Times New Roman','fontsize',16 );
            zlabel(app.UIAxes2,'Z (m)','fontname','Times New Roman','fontsize',16 );
            axis (app.UIAxes2,"equal");
            view(app.UIAxes2,45,10);
            title(app.UIAxes2,'Predicted results')
            StrArrayNew={'Group joint detection result displayed and the points belonging to the same group were assigned to same color！'};
            logRefresh_func(app,StrArrayNew)
        end

        % Callback function
        function SavesamplesButtonPushed(app, event)
             [file,path] = uiputfile({'*.mat','mat-files(*.mat)'},'save');
             str=fullfile(path,file);
             temp=app.pcLearn;
             save(str,'temp');  
        end

        % Menu selected function: InstructionMenu
        function InstructionMenuSelected(app, event)

d = uifigure('Position',[300 300 800 400],'Name','Instruction');
txa = uitextarea(d,...
    'Position',[0 0 800 400],...
    'Value', {'1. Prerequisites for Deployment ';...
    'Verify that version 9.9 (R2020b) of the MATLAB Runtime is installed.   If not, you can run the MATLAB Runtime installer.';...
    'To find its location, enter';...
    '';...
    '    >>mcrinstaller';...
    '';...
    'at the MATLAB prompt.';...
    'NOTE: You will need administrator rights to run the MATLAB Runtime installer. ';...
    '';...
    'Alternatively, download and install the Windows version of the MATLAB Runtime for R2020b from the following link on the MathWorks website:';...
    '';...
    '    https://www.mathworks.com/products/compiler/mcr/index.html';...
    '';...
    'For more information about the MATLAB Runtime and the MATLAB Runtime installer, see "Distribute Applications" in the MATLAB Compiler documentation  in the MathWorks Documentation Center.';...
    '';...
    '2. Files to Deploy and Package';...
    'Files to Package for Standalone ';...
    '================================';...
    '-DisDetANN.exe';...
    '-MCRInstaller.exe ';...
    '    Note: if end users are unable to download the MATLAB Runtime using theinstructions in the previous section, include it when building your component by clicking the "Runtime included in package" link in the Deployment Tool';...
    '';...
    '';...
    '';...
    '3. Definitions';...
    '';...
    'For information on deployment terminology, go to https://www.mathworks.com/help and select MATLAB Compiler >Getting Started > About Application Deployment >Deployment Product Terms in the MathWorks Documentation Center.';...
    '';...
    '';...
    '';...
    '4. Instruction';...
    'Import the raw point cloud data by clicking the File button. It only support.mat/.xls/.xlsx/.pcd/.txt  data .';...
    '';...
    'Parameter section : cited from section 2.3 of the paper';...
    '';...
    '1) Loading';...
    'The filePath will display the file import path.';...
    '';...
    '2) Parameter';...
    '    Click "View raw data" to display the raw point cloud data.';...
    '    Click "XYZ-coordinates" to take the original coordinates as one of the parameters.';...
    '    Type K before clicking "Point normal" and "Point curvature". then click "Point normal" and "Point curvature" respectively to calculate Point normal and Point curvature and take them as one of the parameters.';...
    '';...
    '    Click "Mean distance" to calculate mean distance.';...
    '';...
    '    Type radius  in relation mean distance that before clicking "Point density". radius is generally above or below the mean distance and should not be much greater than mean distance.';...
    '    Click "Point density" to calculate point density  take it as one of the parameters.';...
    '';...
    '3) Plot';...
    '    Which button selected will display point cloud of the color represented by that button value. for example, it  will display point cloud of the color represented by that point normal when user select "Point normal"';...
    '';...
    '';...
    'Network section : cited from section 2.2 and 2.4 of the paper';...
    '';...
    '1) normalization';...
    '    In Network section, the first thing need to do is to click "normalization"  that will normalize the data.';...
    '';...
    '2) Select the samples';...
    '    First, the number of groups on the outcrop must be defined.';...
    '';...
    '    Click "select" to select sample points for each joint group.';...
    '';...
    '    When click "select", need to enter the number of samples user want to select for the joint user have chose in the pop-up dialog box. then, In the pop-up figure, select the points belonging to the same joint through the "data prompt" in the "tool", and click "Enter" on the keyboard for each point selected until the selection is completed.';...
    '';...
    '3) Train the network';...
    '';...
    '    Network include three input arguments that can input from "HiddenSizes"、"PerformFun" and "TrainFun".';...
    '';...
    '    Click "Train" to train the network after input three input arguments.';...
    '';...
    '    Click "Prediction" to predict point cloud data for extracting different sets of joints after the network training is completed.';...
    '    Click "result" to view the point cloud classification results. ';...
    '    If the result is good, then can click "Save net" to save the network, so when open the program again for prediction, user can import the previously saved network by clicking "Load net" without selecting the sample and training the network.';...
    '';...
    '4) Plot frame';...
    '';...
    '    Click "samples" to view the samples user selected.';...
    '';...
    '';...
    'Clustering section : cited from section 2.5 of the paper';...
    '1) Clustering';...
    '';...
    '    DBSCAN algorithm is employed in this program. ';...
    '';...
    '    It includes two parameters, the radius of the neighborhood(ε) and the minimum number of points forming a single joint (minPts), user can input by entering "ε" and "minPts".';...
    '    Click "Cluster" to get individual joint that belonging to the same group joint  user chose in "Group" after input ε and minPts ';...
    '';...
    '    ε is generally above or below the mean distance and should not be much greater than mean distance.';...
    '';...
    '    minPts: user need to observe the clustering results to see if the minPts is appropriate and adjust minPts until the clustering results are satisfactory';...
    '';...
    '2) Plot frame';...
    '';...
    '    Click "Result" to see the clustering results.';...
    '';...
    '';...
    'Orientation section : cited from section 2.6 of the paper';...
    '';...
    '    Click "Calculation" to calculate the orientation of the group joint user chose.';...
    '';...
    '    Click "Result" to see the orientation distribution of the outcrop displayed in "Stereographic Projection".';...
    '';...
    '    If user want to see the orientation value, can click "Export" to export the Excel file.';...
    '';...
    '';...
     '';...
    '';...
    '';...
    '';...
    '';});
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 729 540];
            app.UIFigure.Name = 'MATLAB App';

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.MenuSelectedFcn = createCallbackFcn(app, @FileMenuSelected, true);
            app.FileMenu.Text = '      File     ';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.UIFigure);
            app.HelpMenu.Text = '    Help    ';

            % Create AboutMenu
            app.AboutMenu = uimenu(app.HelpMenu);
            app.AboutMenu.MenuSelectedFcn = createCallbackFcn(app, @AboutMenuSelected, true);
            app.AboutMenu.Text = 'About';

            % Create InstructionMenu
            app.InstructionMenu = uimenu(app.HelpMenu);
            app.InstructionMenu.MenuSelectedFcn = createCallbackFcn(app, @InstructionMenuSelected, true);
            app.InstructionMenu.Text = 'Instruction';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 1 729 540];

            % Create ParameterTab
            app.ParameterTab = uitab(app.TabGroup);
            app.ParameterTab.Title = 'Parameter';
            app.ParameterTab.BackgroundColor = [0.902 0.902 0.902];

            % Create LogPanel
            app.LogPanel = uipanel(app.ParameterTab);
            app.LogPanel.Title = 'Log';
            app.LogPanel.BackgroundColor = [0.902 0.902 0.902];
            app.LogPanel.FontAngle = 'italic';
            app.LogPanel.Position = [11 14 703 106];

            % Create TextArea
            app.TextArea = uitextarea(app.LogPanel);
            app.TextArea.FontSize = 14;
            app.TextArea.Position = [1 0 702 87];

            % Create PlotframePanel
            app.PlotframePanel = uipanel(app.ParameterTab);
            app.PlotframePanel.Title = 'Plot frame';
            app.PlotframePanel.BackgroundColor = [0.902 0.902 0.902];
            app.PlotframePanel.FontAngle = 'italic';
            app.PlotframePanel.FontSize = 13;
            app.PlotframePanel.Position = [347 129 367 381];

            % Create PlotButtonGroup
            app.PlotButtonGroup = uibuttongroup(app.PlotframePanel);
            app.PlotButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @PlotButtonGroupSelectionChanged, true);
            app.PlotButtonGroup.Title = 'Plot';
            app.PlotButtonGroup.FontAngle = 'italic';
            app.PlotButtonGroup.Position = [62 6 252 111];

            % Create PlotXYZcoordinatesButton
            app.PlotXYZcoordinatesButton = uiradiobutton(app.PlotButtonGroup);
            app.PlotXYZcoordinatesButton.Text = 'XYZ-coordinates';
            app.PlotXYZcoordinatesButton.Position = [24 44 111 22];
            app.PlotXYZcoordinatesButton.Value = true;

            % Create PlotPointNormalButton
            app.PlotPointNormalButton = uiradiobutton(app.PlotButtonGroup);
            app.PlotPointNormalButton.Text = 'Point Normal';
            app.PlotPointNormalButton.Position = [139 44 91 22];

            % Create PlotPointCurvatureButton
            app.PlotPointCurvatureButton = uiradiobutton(app.PlotButtonGroup);
            app.PlotPointCurvatureButton.Text = 'Point Curvature';
            app.PlotPointCurvatureButton.Position = [24 15 105 22];

            % Create PlotPointDensityButton
            app.PlotPointDensityButton = uiradiobutton(app.PlotButtonGroup);
            app.PlotPointDensityButton.Text = 'Point Density';
            app.PlotPointDensityButton.Position = [139 15 93 22];

            % Create UIAxes
            app.UIAxes = uiaxes(app.PlotframePanel);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.PlotBoxAspectRatio = [1.33333333333333 1 1];
            app.UIAxes.Position = [6 123 360 229];

            % Create ParameterPanel
            app.ParameterPanel = uipanel(app.ParameterTab);
            app.ParameterPanel.AutoResizeChildren = 'off';
            app.ParameterPanel.Title = 'Parameter ';
            app.ParameterPanel.BackgroundColor = [0.902 0.902 0.902];
            app.ParameterPanel.FontAngle = 'italic';
            app.ParameterPanel.FontSize = 13;
            app.ParameterPanel.Position = [10 129 330 245];

            % Create PointDensityButton
            app.PointDensityButton = uibutton(app.ParameterPanel, 'push');
            app.PointDensityButton.ButtonPushedFcn = createCallbackFcn(app, @PointDensityButtonPushed, true);
            app.PointDensityButton.Position = [24 16 106 22];
            app.PointDensityButton.Text = 'Point Density';

            % Create XYZcoordinatesButton
            app.XYZcoordinatesButton = uibutton(app.ParameterPanel, 'push');
            app.XYZcoordinatesButton.ButtonPushedFcn = createCallbackFcn(app, @XYZcoordinatesButtonPushed, true);
            app.XYZcoordinatesButton.Position = [25 184 106 22];
            app.XYZcoordinatesButton.Text = 'XYZ-coordinates';

            % Create PointnormalButton
            app.PointnormalButton = uibutton(app.ParameterPanel, 'state');
            app.PointnormalButton.ValueChangedFcn = createCallbackFcn(app, @PointnormalButtonValueChanged, true);
            app.PointnormalButton.Text = 'Point normal';
            app.PointnormalButton.Position = [24 142 106 22];

            % Create PointcurvatureButton
            app.PointcurvatureButton = uibutton(app.ParameterPanel, 'state');
            app.PointcurvatureButton.ValueChangedFcn = createCallbackFcn(app, @PointcurvatureButtonValueChanged, true);
            app.PointcurvatureButton.Text = 'Point curvature';
            app.PointcurvatureButton.Position = [24 100 106 22];

            % Create kEditFieldLabel
            app.kEditFieldLabel = uilabel(app.ParameterPanel);
            app.kEditFieldLabel.HorizontalAlignment = 'right';
            app.kEditFieldLabel.Position = [161 142 25 22];
            app.kEditFieldLabel.Text = 'k';

            % Create kEditField
            app.kEditField = uieditfield(app.ParameterPanel, 'numeric');
            app.kEditField.Position = [200 142 60 22];

            % Create radiusEditFieldLabel
            app.radiusEditFieldLabel = uilabel(app.ParameterPanel);
            app.radiusEditFieldLabel.HorizontalAlignment = 'right';
            app.radiusEditFieldLabel.Position = [147 16 38 22];
            app.radiusEditFieldLabel.Text = 'radius';

            % Create radiusEditField
            app.radiusEditField = uieditfield(app.ParameterPanel, 'numeric');
            app.radiusEditField.Position = [200 16 60 22];

            % Create MeandistanceButton
            app.MeandistanceButton = uibutton(app.ParameterPanel, 'push');
            app.MeandistanceButton.ButtonPushedFcn = createCallbackFcn(app, @MeandistanceButtonPushed, true);
            app.MeandistanceButton.Position = [24 58 106 22];
            app.MeandistanceButton.Text = 'Mean distance';

            % Create meandistacneEditField
            app.meandistacneEditField = uieditfield(app.ParameterPanel, 'numeric');
            app.meandistacneEditField.Position = [200 58 60 22];

            % Create LoadingPanel
            app.LoadingPanel = uipanel(app.ParameterTab);
            app.LoadingPanel.Title = 'Loading';
            app.LoadingPanel.BackgroundColor = [0.902 0.902 0.902];
            app.LoadingPanel.FontAngle = 'italic';
            app.LoadingPanel.FontSize = 13;
            app.LoadingPanel.Position = [10 379 330 131];

            % Create ViewrawdataButton
            app.ViewrawdataButton = uibutton(app.LoadingPanel, 'push');
            app.ViewrawdataButton.ButtonPushedFcn = createCallbackFcn(app, @ViewrawdataButtonPushed, true);
            app.ViewrawdataButton.Position = [10 21 91 22];
            app.ViewrawdataButton.Text = 'View raw data';

            % Create FilepathTextAreaLabel
            app.FilepathTextAreaLabel = uilabel(app.LoadingPanel);
            app.FilepathTextAreaLabel.HorizontalAlignment = 'right';
            app.FilepathTextAreaLabel.Position = [15 60 51 22];
            app.FilepathTextAreaLabel.Text = 'File path';

            % Create FilepathTextArea
            app.FilepathTextArea = uitextarea(app.LoadingPanel);
            app.FilepathTextArea.Position = [81 56 244 28];

            % Create NetworkTab
            app.NetworkTab = uitab(app.TabGroup);
            app.NetworkTab.Title = '     Network  ';
            app.NetworkTab.BackgroundColor = [0.902 0.902 0.902];

            % Create logPanel
            app.logPanel = uipanel(app.NetworkTab);
            app.logPanel.Title = 'log';
            app.logPanel.BackgroundColor = [0.902 0.902 0.902];
            app.logPanel.FontAngle = 'italic';
            app.logPanel.FontSize = 13;
            app.logPanel.Position = [6 5 709 103];

            % Create TextArea_2
            app.TextArea_2 = uitextarea(app.logPanel);
            app.TextArea_2.FontSize = 14;
            app.TextArea_2.Position = [1 0 708 84];

            % Create TraintheNetworkPanel
            app.TraintheNetworkPanel = uipanel(app.NetworkTab);
            app.TraintheNetworkPanel.Title = 'Train the Network';
            app.TraintheNetworkPanel.BackgroundColor = [0.902 0.902 0.902];
            app.TraintheNetworkPanel.FontAngle = 'italic';
            app.TraintheNetworkPanel.FontSize = 13;
            app.TraintheNetworkPanel.Position = [7 115 337 193];

            % Create Panel_2
            app.Panel_2 = uipanel(app.TraintheNetworkPanel);
            app.Panel_2.BackgroundColor = [0.902 0.902 0.902];
            app.Panel_2.Position = [0 57 337 116];

            % Create TrainButton
            app.TrainButton = uibutton(app.Panel_2, 'state');
            app.TrainButton.ValueChangedFcn = createCallbackFcn(app, @TrainButtonValueChanged, true);
            app.TrainButton.Text = 'Train';
            app.TrainButton.Position = [219 72 69 22];

            % Create SavenetButton
            app.SavenetButton = uibutton(app.Panel_2, 'push');
            app.SavenetButton.ButtonPushedFcn = createCallbackFcn(app, @SavenetButtonPushed, true);
            app.SavenetButton.Position = [219 30 69 22];
            app.SavenetButton.Text = 'Save net';

            % Create TrainFunDropDownLabel
            app.TrainFunDropDownLabel = uilabel(app.Panel_2);
            app.TrainFunDropDownLabel.HorizontalAlignment = 'right';
            app.TrainFunDropDownLabel.Position = [13 12 53 22];
            app.TrainFunDropDownLabel.Text = 'TrainFun';

            % Create TrainFunDropDown
            app.TrainFunDropDown = uidropdown(app.Panel_2);
            app.TrainFunDropDown.Items = {'trainlm', 'trainbr', 'trainbfg', 'trainrp', 'trainscg', 'traincgb', 'traincgf', 'traincgp', 'trainoss', 'traingdx', 'traingdm', 'traingd'};
            app.TrainFunDropDown.Position = [97 12 69 22];
            app.TrainFunDropDown.Value = 'trainlm';

            % Create PerformFunDropDownLabel
            app.PerformFunDropDownLabel = uilabel(app.Panel_2);
            app.PerformFunDropDownLabel.HorizontalAlignment = 'right';
            app.PerformFunDropDownLabel.Position = [13 47 69 22];
            app.PerformFunDropDownLabel.Text = 'PerformFun';

            % Create PerformFunDropDown
            app.PerformFunDropDown = uidropdown(app.Panel_2);
            app.PerformFunDropDown.Items = {'crossentropy', 'mae', 'mse', 'sae', 'sse', 'msesparse'};
            app.PerformFunDropDown.Position = [97 47 69 22];
            app.PerformFunDropDown.Value = 'crossentropy';

            % Create HiddenSizesEditFieldLabel
            app.HiddenSizesEditFieldLabel = uilabel(app.Panel_2);
            app.HiddenSizesEditFieldLabel.HorizontalAlignment = 'right';
            app.HiddenSizesEditFieldLabel.Position = [13 83 73 22];
            app.HiddenSizesEditFieldLabel.Text = 'HiddenSizes';

            % Create HiddenSizesEditField
            app.HiddenSizesEditField = uieditfield(app.Panel_2, 'numeric');
            app.HiddenSizesEditField.Position = [121 83 43 22];

            % Create Panel_3
            app.Panel_3 = uipanel(app.TraintheNetworkPanel);
            app.Panel_3.BackgroundColor = [0.902 0.902 0.902];
            app.Panel_3.Position = [0 0 178 59];

            % Create LoadnetButton
            app.LoadnetButton = uibutton(app.Panel_3, 'push');
            app.LoadnetButton.ButtonPushedFcn = createCallbackFcn(app, @LoadnetButtonPushed, true);
            app.LoadnetButton.BackgroundColor = [0.9412 0.9412 0.9412];
            app.LoadnetButton.Position = [37 19 100 22];
            app.LoadnetButton.Text = 'Load net';

            % Create Panel_4
            app.Panel_4 = uipanel(app.TraintheNetworkPanel);
            app.Panel_4.BackgroundColor = [0.902 0.902 0.902];
            app.Panel_4.Position = [176 0 161 59];

            % Create PredictionButton
            app.PredictionButton = uibutton(app.Panel_4, 'state');
            app.PredictionButton.ValueChangedFcn = createCallbackFcn(app, @PredictionButtonValueChanged, true);
            app.PredictionButton.Text = 'Prediction';
            app.PredictionButton.Position = [44 18 69 22];

            % Create SelecttheSamplesPanel
            app.SelecttheSamplesPanel = uipanel(app.NetworkTab);
            app.SelecttheSamplesPanel.Title = 'Select the Samples ';
            app.SelecttheSamplesPanel.BackgroundColor = [0.902 0.902 0.902];
            app.SelecttheSamplesPanel.FontAngle = 'italic';
            app.SelecttheSamplesPanel.FontSize = 13;
            app.SelecttheSamplesPanel.Position = [7 315 337 107];

            % Create GroupsinthisoutcropEditFieldLabel
            app.GroupsinthisoutcropEditFieldLabel = uilabel(app.SelecttheSamplesPanel);
            app.GroupsinthisoutcropEditFieldLabel.HorizontalAlignment = 'right';
            app.GroupsinthisoutcropEditFieldLabel.Position = [9 47 126 22];
            app.GroupsinthisoutcropEditFieldLabel.Text = ' Groups in this outcrop';

            % Create GroupsinthisoutcropEditField
            app.GroupsinthisoutcropEditField = uieditfield(app.SelecttheSamplesPanel, 'numeric');
            app.GroupsinthisoutcropEditField.ValueChangedFcn = createCallbackFcn(app, @GroupsinthisoutcropEditFieldValueChanged, true);
            app.GroupsinthisoutcropEditField.Position = [231 47 60 22];

            % Create GroupDropDownLabel
            app.GroupDropDownLabel = uilabel(app.SelecttheSamplesPanel);
            app.GroupDropDownLabel.HorizontalAlignment = 'right';
            app.GroupDropDownLabel.Position = [13 15 39 22];
            app.GroupDropDownLabel.Text = 'Group';

            % Create GroupDropDown
            app.GroupDropDown = uidropdown(app.SelecttheSamplesPanel);
            app.GroupDropDown.Items = {'Option 1', 'Option 2'};
            app.GroupDropDown.Position = [67 15 100 22];

            % Create SelectButton
            app.SelectButton = uibutton(app.SelecttheSamplesPanel, 'push');
            app.SelectButton.ButtonPushedFcn = createCallbackFcn(app, @SelectButtonPushed, true);
            app.SelectButton.Position = [231 15 60 22];
            app.SelectButton.Text = 'Select';

            % Create PlotframePanel_2
            app.PlotframePanel_2 = uipanel(app.NetworkTab);
            app.PlotframePanel_2.Title = 'Plot frame';
            app.PlotframePanel_2.BackgroundColor = [0.902 0.902 0.902];
            app.PlotframePanel_2.FontAngle = 'italic';
            app.PlotframePanel_2.FontSize = 13;
            app.PlotframePanel_2.Position = [354 115 361 389];

            % Create SamplesButton_2
            app.SamplesButton_2 = uibutton(app.PlotframePanel_2, 'push');
            app.SamplesButton_2.ButtonPushedFcn = createCallbackFcn(app, @SamplesButton_2Pushed, true);
            app.SamplesButton_2.Position = [56 48 100 22];
            app.SamplesButton_2.Text = 'Samples';

            % Create ResultButton_4
            app.ResultButton_4 = uibutton(app.PlotframePanel_2, 'push');
            app.ResultButton_4.ButtonPushedFcn = createCallbackFcn(app, @ResultButton_4Pushed, true);
            app.ResultButton_4.Position = [207 48 100 22];
            app.ResultButton_4.Text = 'Result';

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.PlotframePanel_2);
            title(app.UIAxes2, 'Title')
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.PlotBoxAspectRatio = [1.47058823529412 1 1];
            app.UIAxes2.Position = [27 113 324 247];

            % Create NormalizationPanel
            app.NormalizationPanel = uipanel(app.NetworkTab);
            app.NormalizationPanel.Title = 'Normalization';
            app.NormalizationPanel.BackgroundColor = [0.902 0.902 0.902];
            app.NormalizationPanel.FontAngle = 'italic';
            app.NormalizationPanel.FontSize = 13;
            app.NormalizationPanel.Position = [7 433 337 71];

            % Create normaliztionButton
            app.normaliztionButton = uibutton(app.NormalizationPanel, 'push');
            app.normaliztionButton.ButtonPushedFcn = createCallbackFcn(app, @normaliztionButtonPushed, true);
            app.normaliztionButton.Position = [115 15 100 22];
            app.normaliztionButton.Text = 'normaliztion';

            % Create ClusteringTab
            app.ClusteringTab = uitab(app.TabGroup);
            app.ClusteringTab.Title = '   Clustering';
            app.ClusteringTab.BackgroundColor = [0.902 0.902 0.902];

            % Create ClusteringPanel
            app.ClusteringPanel = uipanel(app.ClusteringTab);
            app.ClusteringPanel.Title = 'Clustering';
            app.ClusteringPanel.BackgroundColor = [0.902 0.902 0.902];
            app.ClusteringPanel.FontAngle = 'italic';
            app.ClusteringPanel.FontSize = 13;
            app.ClusteringPanel.Position = [10 145 302 358];

            % Create GroupDropDown_2Label
            app.GroupDropDown_2Label = uilabel(app.ClusteringPanel);
            app.GroupDropDown_2Label.HorizontalAlignment = 'right';
            app.GroupDropDown_2Label.Position = [9 153 39 22];
            app.GroupDropDown_2Label.Text = 'Group';

            % Create GroupDropDown_2
            app.GroupDropDown_2 = uidropdown(app.ClusteringPanel);
            app.GroupDropDown_2.Items = {'Option 1', 'Option 2'};
            app.GroupDropDown_2.Position = [107 153 93 22];

            % Create minPtsEditFieldLabel
            app.minPtsEditFieldLabel = uilabel(app.ClusteringPanel);
            app.minPtsEditFieldLabel.HorizontalAlignment = 'right';
            app.minPtsEditFieldLabel.Position = [9 191 42 22];
            app.minPtsEditFieldLabel.Text = 'minPts';

            % Create minPtsEditField
            app.minPtsEditField = uieditfield(app.ClusteringPanel, 'numeric');
            app.minPtsEditField.Position = [136 191 63 22];

            % Create Label
            app.Label = uilabel(app.ClusteringPanel);
            app.Label.HorizontalAlignment = 'right';
            app.Label.FontSize = 13;
            app.Label.Position = [9 228 25 22];
            app.Label.Text = 'ε';

            % Create epsEditField
            app.epsEditField = uieditfield(app.ClusteringPanel, 'numeric');
            app.epsEditField.Position = [136 228 63 22];

            % Create ClusterButton
            app.ClusterButton = uibutton(app.ClusteringPanel, 'push');
            app.ClusterButton.ButtonPushedFcn = createCallbackFcn(app, @ClusterButtonPushed, true);
            app.ClusterButton.Position = [224 153 64 22];
            app.ClusterButton.Text = 'Cluster';

            % Create ResultButton_2
            app.ResultButton_2 = uibutton(app.ClusteringPanel, 'push');
            app.ResultButton_2.ButtonPushedFcn = createCallbackFcn(app, @ResultButton_2Pushed, true);
            app.ResultButton_2.Position = [107 110 94 22];
            app.ResultButton_2.Text = 'Result ';

            % Create LogPanel_2
            app.LogPanel_2 = uipanel(app.ClusteringTab);
            app.LogPanel_2.Title = 'Log';
            app.LogPanel_2.BackgroundColor = [0.902 0.902 0.902];
            app.LogPanel_2.FontAngle = 'italic';
            app.LogPanel_2.FontSize = 13;
            app.LogPanel_2.Position = [10 13 704 113];

            % Create TextArea_3
            app.TextArea_3 = uitextarea(app.LogPanel_2);
            app.TextArea_3.FontSize = 14;
            app.TextArea_3.Position = [2 0 702 91];

            % Create PotframePanel
            app.PotframePanel = uipanel(app.ClusteringTab);
            app.PotframePanel.Title = 'Pot frame';
            app.PotframePanel.BackgroundColor = [0.902 0.902 0.902];
            app.PotframePanel.FontAngle = 'italic';
            app.PotframePanel.FontSize = 13;
            app.PotframePanel.Position = [325 145 390 358];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.PotframePanel);
            title(app.UIAxes3, 'Title')
            xlabel(app.UIAxes3, 'X')
            ylabel(app.UIAxes3, 'Y')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.PlotBoxAspectRatio = [1.33064516129032 1 1];
            app.UIAxes3.Position = [22 18 348 289];

            % Create OrientationTab
            app.OrientationTab = uitab(app.TabGroup);
            app.OrientationTab.Title = '  Orientation';
            app.OrientationTab.BackgroundColor = [0.902 0.902 0.902];

            % Create ExportButton
            app.ExportButton = uibutton(app.OrientationTab, 'push');
            app.ExportButton.ButtonPushedFcn = createCallbackFcn(app, @ExportButtonPushed, true);
            app.ExportButton.Position = [64 226 100 22];
            app.ExportButton.Text = 'Export';

            % Create ResultButton_3
            app.ResultButton_3 = uibutton(app.OrientationTab, 'push');
            app.ResultButton_3.ButtonPushedFcn = createCallbackFcn(app, @ResultButton_3Pushed, true);
            app.ResultButton_3.Position = [64 287 100 22];
            app.ResultButton_3.Text = 'Result';

            % Create CalculationButton
            app.CalculationButton = uibutton(app.OrientationTab, 'push');
            app.CalculationButton.ButtonPushedFcn = createCallbackFcn(app, @CalculationButtonPushed, true);
            app.CalculationButton.Position = [195 359 100 22];
            app.CalculationButton.Text = 'Calculation';

            % Create StereographicProjectionPanel
            app.StereographicProjectionPanel = uipanel(app.OrientationTab);
            app.StereographicProjectionPanel.Title = 'Stereographic Projection';
            app.StereographicProjectionPanel.BackgroundColor = [0.902 0.902 0.902];
            app.StereographicProjectionPanel.FontAngle = 'italic';
            app.StereographicProjectionPanel.FontSize = 13;
            app.StereographicProjectionPanel.Position = [334 129 379 374];

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.StereographicProjectionPanel);
            title(app.UIAxes4, 'Title')
            xlabel(app.UIAxes4, 'X')
            ylabel(app.UIAxes4, 'Y')
            zlabel(app.UIAxes4, 'Z')
            app.UIAxes4.PlotBoxAspectRatio = [1.22393822393822 1 1];
            app.UIAxes4.Position = [19 24 351 290];

            % Create LogPanel_3
            app.LogPanel_3 = uipanel(app.OrientationTab);
            app.LogPanel_3.Title = 'Log';
            app.LogPanel_3.FontAngle = 'italic';
            app.LogPanel_3.FontSize = 13;
            app.LogPanel_3.Position = [10 17 703 103];

            % Create TextArea_4
            app.TextArea_4 = uitextarea(app.LogPanel_3);
            app.TextArea_4.FontSize = 14;
            app.TextArea_4.Position = [1 -2 702 83];

            % Create GroupDropDown_3
            app.GroupDropDown_3 = uidropdown(app.OrientationTab);
            app.GroupDropDown_3.Items = {'Option 1', 'Option 2'};
            app.GroupDropDown_3.Position = [64 359 100 22];

            % Create GroupDropDown_3Label
            app.GroupDropDown_3Label = uilabel(app.OrientationTab);
            app.GroupDropDown_3Label.HorizontalAlignment = 'right';
            app.GroupDropDown_3Label.Position = [10 359 39 22];
            app.GroupDropDown_3Label.Text = 'Group';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DisDetANN

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end