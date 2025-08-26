 function dataOut = getErrors(xdata,ydata)
        %function dataOut = getErrors(xdata,ydata)
        %From this web site:
        %http://terpconnect.umd.edu/~toh/spectrum/LeastSquaresMatlab.txt
        %KL modifying 4/21/2014
        
        x = xdata;
        y = ydata;
        % Simple Matlab script for calculating the first-order least-square fit of y vs x,
        % including the Slope and Intercept and the predicted standard deviation of
        % the slope (SDSlope) and intercept (SDIntercept).
        
        NumPoints=length(x);
        Sxx = sum((x-mean(x)).^2);
        Syy = sum((y-mean(y)).^2);
        Sxy = sum((x-mean(x)).*(y-mean(y)));
        Slope = Sxy./Sxx;
        Intercept = mean(y)-Slope*mean(x);
        
        Sy = sqrt((Syy-Slope^2*Sxx)/(NumPoints-2));
        
        SDslope = Sy/sqrt(Sxx);
        SDintercept = Sy*sqrt(1./(NumPoints-(sum(x).^2)./sum(x.^2)));
        
        r2 = 1 - ((Syy-Slope^2*Sxx) ./Syy);
        
        %data to send out of this function (when it is a function)
        dataOut.slope = Slope;
        dataOut.intercept = Intercept;
        dataOut.SDslope = SDslope;
        dataOut.PercentSlopeError = SDslope./Slope;
        dataOut.SDintercept = SDintercept;
        dataOut.PercentInterceptError = SDintercept./Intercept;
        dataOut.r2 = r2;
        dataOut.Sxx = Sxx;
        dataOut.Syy = Syy;
        dataOut.Sxy = Sxy;
        dataOut.Sy = Sy;
        
        % Now to calculate prediction intervals at 95% confidence.
        A = (1-r2)*((NumPoints-1)/(NumPoints-2))*(1/Sxx);
        B = tinv(0.975,NumPoints-2)*Sy;
        C = Sxx*(1+(1/NumPoints));
        xM = mean(xdata);
        dataOut.A = A;
        dataOut.B = B;
        dataOut.C = C;
        dataOut.xM = xM;
        
    end %end of getErrors as a function