%Lidar sensor which can be used to extract measurements from trajectories,
%clutter can be added as well.
classdef LidarSensor
    properties
        maxRange
        noiseStd
        bearingAngles
        clutterIntensity
        minRangeClutter
    end

    methods
        function obj = LidarSensor(maxRange, noiseStd, dbearing, clutterIntensity, minRangeClutter)
            obj.maxRange = maxRange;
            obj.noiseStd = noiseStd;
            obj.bearingAngles = 0:dbearing:2*pi-dbearing;
            obj.clutterIntensity = clutterIntensity;
            obj.minRangeClutter = minRangeClutter;
        end
        
        function positions = getMeasurements(obj, shipTrajectories, addClutter)
            Nbearings = length(obj.bearingAngles);
            Nship = length(shipTrajectories);
            Nt = size(shipTrajectories{1}.xKin,2); % TODO: Assumes Nship>=1       
            positions = cell(Nt,1);
            ranges = cell(Nt,1);
            bearings = cell(Nt,1);
%             bearing_indexes = cell(Nt,1);            
            for nt = 1:Nt
                rangeComparison = inf(Nbearings, Nship+1);
                for ns=1:Nship
                    [extendRanges, extendBearingIndexes] = ...
                        obj.intersectRaysWithClosedCurve(shipTrajectories{ns}.extendNE(:,:,nt));
                    extendRanges = extendRanges + obj.noiseStd*randn(length(extendRanges),1);
                    rangeComparison(extendBearingIndexes, ns) = extendRanges;
                end
                if addClutter
                    Nclutter = round(obj.clutterIntensity*pi*obj.maxRange^2);
                    clutterRanges = (obj.maxRange-obj.minRangeClutter)*rand(min(Nclutter,Nbearings),1)+obj.minRangeClutter;
                    clutterBearingIndexes = randperm(Nbearings, min(Nclutter,Nbearings));
                    rangeComparison(clutterBearingIndexes, end) = clutterRanges;
                end
%                 if addPhantomShip
%                     max_rotation = pi/2;
%                     max_measurements = 50;
%                     psi = max_rotation*rand(1);
%                     hull_center = mean(hulls_n(:,:,nt),2);
%                     phantom_hull = RotAndTrans(hulls_n(:,:,nt)-hull_center, psi, hull_center);
%                     [phantom_ranges, phantom_bearing_indexes] = ...
%                     intersectRaysClosedCurve(obj.bearings, phantom_hull);
%                     Nphantom = min(max_measurements, length(phantom_ranges));
%                     Nphantom = randi(Nphantom,1);
%                     selected_indexes = randperm(length(phantom_ranges), Nphantom);
%                     phantom_bearing_indexes = phantom_bearing_indexes(selected_indexes);
%                     range_comparison(phantom_bearing_indexes, 3) = phantom_ranges(selected_indexes);
%                 end
                [rangeComparisonMin, ~] = min(rangeComparison, [], 2);
                measurementIndexes = ~isinf(rangeComparisonMin);
                ranges{nt} = rangeComparisonMin(measurementIndexes);
                bearings{nt} = obj.bearingAngles(measurementIndexes)';
                [N, E] = pol2cart(bearings{nt}, ranges{nt});
                positions{nt} = [N'; E'];
%                 bearing_indexes{nt} = find(measurement_indexes);
            end
        end
        
        function [ranges, bearingIndexes] = intersectRaysWithClosedCurve(obj, closedCurve)
            Nbearing = length(obj.bearingAngles);
            ranges = inf(Nbearing, 1);
            Npoint = size(closedCurve, 2);
            for nb = 1:Nbearing
                bearing = obj.bearingAngles(nb);
                dirRay = [cos(bearing); sin(bearing)];
                dirRayPerp = [-sin(bearing); cos(bearing)];
                for np = 1:Npoint
                    x1 = closedCurve(1,np);
                    y1 = closedCurve(2,np);
                    if np < Npoint
                        x2 = closedCurve(1,np+1);
                        y2 = closedCurve(2,np+1);
                    else
                        x2 = closedCurve(1,1);
                        y2 = closedCurve(2,1);
                    end
                    if ([x1 y1]*dirRayPerp)*([x2 y2]*dirRayPerp)<=0
                       A = [dirRayPerp(1) dirRayPerp(2);
                           y2-y1 -(x2-x1)];
                       b = [0; (y2-y1)*x1-(x2-x1)*y1];
                       posCross = A\b;
                       if dirRay.'*posCross>=0 && norm(posCross) < ranges(nb)
                           ranges(nb) = norm(posCross);
                       end
                    end
                end
            end
            bearingIndexes = 1:Nbearing;
            criteria = ~isinf(ranges) + ranges<obj.maxRange;
            bearingIndexes = bearingIndexes(criteria)';
            ranges = ranges(criteria);
        end
        function [ranges,bearings] = intersectRays(obj, closedCurve)
            [ranges,bearingIndices] = obj.intersectRaysWithClosedCurve(closedCurve);
            bearings = obj.bearingAngles(bearingIndices);
        end
    end
end