f = uifigure("Name","Linear Algebraic Equations Solver", "NumberTitle","off", "Color","#1772CD");
f.InnerPosition = [200 200 800 600];
f.Resize = "off";
%grid = uigridlayout(f,[4 4]);
format long;
global method;
global stop_c_select;
global message;
global p_info;
global roots;
message = "";
stop_c_select = "";
method = "";
p_info = 0;
roots = [];

%specify the group for the buttons that dictate which method will be used
method_group = uibuttongroup(f, "Position",[0 500 800 100], ...
    "BackgroundColor","#ffffff");
method_group.SelectionChangedFcn = @(method_group, event) set_method(event);
%method_group.SelectionChangedFcn = (@set_method);

%invisible placeholder button so thatthe method buttons can be pressed with
%no errors
p0 = uitogglebutton(method_group, "Text", "");
p0.InnerPosition = [0,0,0,0];

p1 = uitogglebutton(method_group,"Text",'Gaussian Elimnation');
%resize the text box
p1.InnerPosition = [30 10 150 70];
p1.BackgroundColor = "#70CEAB";
p1.FontWeight = "Bold";
p1.FontName = "Arial";
p1.HorizontalAlignment = "Center";
p1.FontSize = 14;

p2 = uitogglebutton(method_group,"Text",'Gauss-Seidel');
p2.InnerPosition = [220 10 150 70];
p2.BackgroundColor = "#70CEAB";
p2.FontWeight = "Bold";
p2.FontName = "Arial";
p2.HorizontalAlignment = "Center";
p2.FontSize = 14;

p3 = uitogglebutton(method_group,"Text",'Jacobi');
p3.InnerPosition = [420 10 150 70];
p3.BackgroundColor = "#70CEAB";
p3.FontWeight = "Bold";
p3.FontName = "Arial";
p3.HorizontalAlignment = "Center";
p3.FontSize = 14;

p4 = uitogglebutton(method_group,"Text",'Cramer');
p4.Position = [440 370 100 40];
p4.InnerPosition = [620 10 150 70];
p4.BackgroundColor = "#70CEAB";
p4.FontWeight = "Bold";
p4.FontName = "Arial";
p4.HorizontalAlignment = "Center";
p4.FontSize = 14;

% n_label = uilabel(f, "Position", [250 450 300 40], "Text", "Enter the size n of the augmented matrix");
% n_label.FontSize = 16;
% n_label.FontColor = "#ffffff";
% 
% n_textbox = uitextarea(f,"Placeholder","Eg: 3","Position",[370 410 50 40]);
% n_textbox.HorizontalAlignment = "center";
% n_textbox.ValueChangedFcn = @(n_textbox,event) set_n(n_textbox);

divider = uilabel(f,"Position",[530 0 5 500] ...
    ,"BackgroundColor","#000000");

left_label = uilabel(f, "Position", [10 450 200 50], "Text", "Enter the values for the left hand side of the equations seperated by commas on different lines");
left_label.FontColor = "#ffffff";
left_label.FontWeight = "Bold";
left_label.WordWrap = "on";
left_label.BackgroundColor = "#000000";


left_textbox = uitextarea(f,"Position",[10 300 200 150], "Placeholder",sprintf("x0,x1,x2\nx3,x4,x5\nx6,x7,x9"));
left_textbox.ValueChangedFcn = @(left_textbox,event) get_matrix_left(left_textbox);
left_textbox.FontSize = 20;
left_textbox.FontWeight = "Bold";

right_label = uilabel(f, "Position", [230 450 70 50], "Text", "Enter values for right hand side");
right_label.FontColor = "#ffffff";
right_label.FontWeight = "Bold";
right_label.WordWrap = "on";
right_label.BackgroundColor = "#000000";

right_textbox = uitextarea(f,"Position",[230 300 70 150], "Placeholder",sprintf("b0\nb1\nb2"));
right_textbox.ValueChangedFcn = @(right_textbox,event) get_matrix_right(right_textbox);
right_textbox.FontSize = 20;
right_textbox.FontWeight = "Bold";

stop_criteria = "";
stop_c_panel = uipanel(f,"Position",[320 300 110 150], "Title","Stopping Criteria", ...
    "FontWeight","bold","BackgroundColor","#21C8B3", "TitlePosition","centertop");
mae_button = uibutton(stop_c_panel, "push","InnerPosition",[20 70 70 50], ...
    "Text","MAE", "FontColor", "#ffffff", "BackgroundColor", "#8A27EB", "FontWeight", "Bold");
mae_button.ButtonPushedFcn = @(mae_button,event) s_criteria(mae_button.Text);

rmse_button = uibutton(stop_c_panel, "push","InnerPosition",[20 10 70 50], ...
    "Text","RMSE", "FontColor", "#ffffff", "BackgroundColor", "#EB9427", "FontWeight", "Bold");
rmse_button.ButtonPushedFcn = @(rsme_button,event) s_criteria(rsme_button.Text);


th_label = uilabel(f, "Position", [440 240 70 50], "Text", "Enter Tolerance Threshold");
th_label.FontColor = "#ffffff";
th_label.FontWeight = "Bold";
th_label.WordWrap = "on";
th_label.BackgroundColor = "#000";

th_textbox = uitextarea(f,"Position",[440 190 70 50], "Placeholder","0.0001");
th_textbox.ValueChangedFcn = @(th_textbox,event) get_threshold(th_textbox);
th_textbox.FontSize = 12;
th_textbox.FontWeight = "Bold";



s_approx_label = uilabel(f, "Position", [440 450 85 50], "Text", "Enter starting approximation");
s_approx_label.FontColor = "#ffffff";
s_approx_label.FontWeight = "Bold";
s_approx_label.WordWrap = "on";
s_approx_label.BackgroundColor = "#000";

s_approx_textbox = uitextarea(f,"Position",[440 300 70 150],"Placeholder", sprintf("s0\ns1\ns2"));
s_approx_textbox.FontSize = 18;
s_approx_textbox.FontWeight = "Bold";

s_approx_textbox.ValueChangedFcn = @(s_approx_textbox,event) get_s_approx(s_approx_textbox);




error_box_label = uilabel(f,"Text", "MESSAGES", ...
    "Position",[0 150 530 30]);
error_box_label.HorizontalAlignment = "Center";
error_box_label.FontWeight = "Bold";
error_box_label.BackgroundColor = "#BFBB13";

error_box = uilabel(f, "Position", [0 0 530 150]);
error_box.BackgroundColor = "#DE6363";
%error_box.Editable = "off";
%error_box.Value = "Hello ";
error_box.FontColor = "#ffffff";
error_box.FontWeight = "Bold";
error_box.FontSize = 14;
error_box.Text = "";
error_box.FontColor = "#ffffff";

%error_box.Enable = "on";




results_panel = uipanel(f, "Position",[540 20 250 400], ...
    "BackgroundColor","#7F6FE1");

roots_results_label = uilabel(results_panel,"Text","ROOTS", ...
    "InnerPosition",[0 370 250 30], ...
    "FontWeight","Bold", "BackgroundColor","#ffffff", ...
    "HorizontalAlignment","Center", "FontSize",18, "FontColor","#44BA1B");
roots_results_box = uilabel(results_panel, ...
    "InnerPosition",[0 220 250 150], "BackgroundColor","#6FE1E0" ...
    , "HorizontalAlignment","Center" ...
    , "FontColor", "#000000", "FontWeight","Bold", "FontSize", 16, "Text", "");

%roots_results_box.ValueChangedFcn = @(main_button,event) display_results(roots_results_box);
roots_results_box.Text = "";


mean_t_err_results_label = uilabel(results_panel,"Text","MEAN TRUE ERROR", ...
    "InnerPosition",[0 180 250 30], ...
    "FontWeight","Bold", "BackgroundColor","#ffffff", ...
    "HorizontalAlignment","Center", "FontSize",12, "FontColor","#44BA1B");
mean_t_err_box = uilabel(results_panel, ...
    "InnerPosition",[0 132 250 50], "BackgroundColor","#175EAB" ...
    , "HorizontalAlignment","Center" ...
    , "FontColor", "#fff", "FontWeight","Bold", "FontSize", 16, "Text", "");

main_button = uibutton(results_panel,"Text", "CALCULATE", ...
    "InnerPosition",[75 0 100 80], "FontWeight","Bold", "FontSize",15 ...
    , "BackgroundColor","#83EF68");

main_button.ButtonPushedFcn = @(main_button,event) check_input(get_matrix_left(left_textbox), get_matrix_right(right_textbox), ...
    method, stop_c_select,get_threshold(th_textbox), get_s_approx(s_approx_textbox),roots_results_box,mean_t_err_box, error_box);



% example_button = uibutton(f, "Position", [700 430 70 50] ...
%     ,"Text","Example");
% example_button.ButtonPushedFcn = @(example_button,event) example(roots_results_box);


% function example(roots_results_box)
%     
%     A = [3 1 -4 ; -2 3 1 ; 2 0 5];
%     b = [7;-5;10];
%     global method;
%     method = "Gauss-Seidel";
%     global stop_c_select;
%     stop_c_select = "MAE";
%     check_input(A, b,method,stop_c_select, 0.0001,[0;0;0], roots_results_box);
% end


function [A] = get_matrix_left(left_textbox)

    arr = left_textbox.Value;

%     if length(arr) ~= n 
%         fprintf("\n\nValue of n is not equal to array")
%     end

    for i = 1 : length(arr)

        %split the entries into string arrays with a comma as the delimeter
        row  = strsplit(string(arr(i)), ",");

        %trim white spaces from all the elements in the row
        row = strtrim(row);
        for j = 1:length(row)
            %convert each string element in the row to a number and assign
            %it a respective place in the matrix
            A(i,j) = str2double(row(j));
            %fprintf("\n %d %d", i, j)
        end
        
    end

end

function [b] = get_matrix_right(right_textbox)

    arr = right_textbox.Value;

    %     if length(arr) ~= n 
    %         fprintf("\n\nValue of n is not equal to array")
    %     end

    for i = 1 : length(arr)

        %split the entries into string arrays with a comma as the delimeter
        row  = strsplit(string(arr(i)));

        %trim white spaces from all the elements in the row
        row = strtrim(row);
        for j = 1:length(row)
            %convert each string element in the row to a number and assign
            %it a respective place in the matrix
            b(i,j) = str2double(row(j));
            %fprintf("\n %d %d", i, j)
        end
        
    end

    %j = b

end

function check_input(A,b,method, stop_c_select, threshold, s_approx,roots_results_box,mean_t_err_box,error_box)
    clear_message();
    roots_results_box.Text = "";
    error_box.Text
    global message;
    global method;
    global stop_c_select;
    global roots;
    check = true;
    check2 = true;



    if size(A,1) ~= size(b,1)
        add_to_message("**ROWS OF MATRIX A & b MUST HAVE THE SAME SIZE***");
        check = false;
        check2 = false;
    end

    if size(A,1) ~= size(A,2)
        add_to_message("**ROWS OF MATRIX A MUST BE OF SIZE nxn***");
        check = false;
        check2 = false;
    end

    if check2 == true
        %check for starting approximation
        if method == "Jacobi" || method == "Gauss-Seidel"
            
            if method == "Jacobi"
               add_to_message("METHOD: " + method);
            end
            
            if method == "Gauss-Seidel"
                add_to_message("METHOD: " + method);
            end
            
            
            if ~isnan(s_approx)
                    add_to_message("STARTING APPROXIMATION: "+mat2str(s_approx));
            end
            
            if ~isnan(threshold)
                add_to_message("THRESHOLD: "+string(threshold));
            end
    
            if isnan(s_approx)
                add_to_message("NO START APPROXIMATION GIVEN. SET TO DEFAULT [0 0 0]");
                s_approx = [0;0;0];
            end
            
            if isnan(threshold)
                add_to_message("NO TOLERANCE THRESHOLD GIVEN. SET TO DEFAULT 0.0001");
                threshold = 0.0001;
            end
            
            
             if check_matrix([A,b]) == 0
                 add_to_message("**NOT DIAGONALLY DOMINANT. A SOLUTION MIGHT NOT BE GENERATED**");
             end
            
             if check_matrix([A,b]) == 1
                 add_to_message("Matrix is diagonally dominant");
             end
            
        end
        
    
       if det(A) == 0
           add_to_message("**The matrix is singular. No solutions can be found**");
           check = false;
       else
           add_to_message("The matrix is not singular");
       end
    end


    if method == "Jacobi" || method == "Gauss-Seidel"
        if isnan(threshold)
            add_to_message("**SPECIFY A THRESHOLD***");
            check = false;
        end
        if stop_c_select == ""
            add_to_message("**SPECIFY A STOPPING CRITERION**");
            check = false;
        end
    end
    
    if method == ""
       add_to_message("**SPECIFY A METHOD**"); 
       check = false;
    end
   
    if check == true && check2 == true
        aug = [A,b];
        
        if method ~= "Jacobi" && method ~= "Gauss-Seidel"
            add_to_message("METHOD: " + method);
        end
        if method == "Jacobi" || method == "Gauss-Seidel"

            add_to_message("STOPPING CRITERION: "+stop_c_select);      
        end
        if method == "Gauss-Seidel"
            res = gauss_seidel(aug,s_approx,threshold,stop_c_select);
        elseif method == "Cramer"
            res = Cramer(aug);
        elseif method == "Jacobi"
            res = jacobi(aug,s_approx,threshold,stop_c_select);
        elseif method == "Gaussian Elimnation"
            res = gaussian(A,b);
        end
        disp(res)%%% 
        disp(isinf(res))%%%
        %try to make diagonally dominant matrix in the reult is -inf
        inf_check = false;
        for i = 1:size(res,1)
            if isinf(res(i)) || isnan(res(i))
                inf_check = true;
            end
        end
        if method == "Jacobi" || method == "Gauss-Seidel"
            if inf_check == true
                add_to_message("<Attempted to make matrix diagonally dominant>*");
                aug2 = make_dd(aug);
                if method == "Jacobi"
                    res = jacobi(aug2,s_approx,threshold,stop_c_select);
                elseif method == "Gauss-Seidel"
                    res = gauss_seidel(aug2,s_approx,threshold,stop_c_select);
                end
            end
        end

        %check for incorrect values attained from incorrect input
        bad_input_check = false;
        for i = 1:size(res,1)
            if isinf(res(i)) || isnan(res(i))
                bad_input_check = true;
            end
        end

        if bad_input_check == true
            add_to_message("**Please make sure that the input is correct**");
        else
            roots = res;
            roots_results_box.Text = string(res);
            mean_t_err_box.Text = string(mean_true_error(A,b,res,length(res)));
        end

    end
    
   
    
   %disp(message);
   error_box.Text = message;

end


function [answers] = gauss_seidel(matrix, guess, tolerance, stopping) %Define the function
    %Algorithm for Gauss-Seidel
    [rows,cols] = size(matrix);
    roots = matrix(:, cols);
    matrix(:, cols) = [];
    [newRows, newCols] = size(matrix);
    
    for i=1:newRows
        roots(i) = roots(i)/matrix(i,i);
        xValues(i) = guess(i);
        oldXValues(i) = xValues(i);
        for j = 1:newRows
            if(i~=j)
                matrix(i,j) = matrix(i,j) / matrix(i,i);
            end
        end
    end
    error = intmax;
    
    iterations = 0;
    while(error > tolerance)
        iterations = iterations + 1;
        error = 0;
        for i = 1:newRows
            xValues(i) = roots(i);
            for j = 1:newRows
                if(i~=j)
                    xValues(i) = xValues(i) - matrix(i,j) * xValues(j);
                end
            end
            error = error + abs(xValues(i) - oldXValues(i));
            oldXValues(i) = xValues(i);
        end
        if stopping == "MAE"
            error = error/newRows;
        elseif stopping == "RMSE"
            error = sqrt((error^2)/newRows);
        end
    end
    answers = xValues;
    i = iterations;
end

function [x] = Cramer(Aug)
    % Using Cramer's rule for solving a system of n equations in n unknowns. 
    % input:
    %   Aug = augmented matrix of the system
    % output:
    %   x = solution vector
    [n,m] = size(Aug);
    if n+1~=m, error('Incorrect Matrix Aug size.'); end
    % check whether the determinant of the system matrix is non-zero
    A = Aug(:,1:n);
    D = det(A);
    x = zeros(n,1);
    if ( D == 0 )
        error('The system is singular. There is no solution.'); 
    else
        for i=1:n
            Ai = A;
            Ai(:,i) = Aug(:,n+1);
            x(i) = det(Ai)/D;
        end
    end
end

function [answers] = jacobi(matrix, guess, tolerance, stopping) %Define the function
    %Algorithm for Jacobi
    [rows,cols] = size(matrix);
    roots = matrix(:, cols);
    matrix(:, cols) = [];
    [newRows, newCols] = size(matrix);
    
    for i = 1:newRows
        roots(i) = roots(i)/matrix(i,i);
        newXValues(i) = guess(i);
        for j = 1:newRows
            if(i~=j)
                matrix(i,j) = matrix(i,j)/matrix(i,i);
            end
        end
    end
    
    error = intmax;
    iterations = 0;
    while(error > tolerance)
        iterations = iterations + 1;
        error = 0;
        for i = 1: newRows
            oldXValue(i) = newXValues(i);
            newXValues(i) = roots(i);
        end
        for i = 1: newRows
            for j = 1: newRows
                if(i~=j)
                    newXValues(i) = newXValues(i)- (matrix(i,j) * oldXValue(j));
                end
            end
            error = error + abs(newXValues(i)-oldXValue(i));
        end
        if stopping == "MAE"
            error = error/newRows;
        elseif stopping == "RMSE"
            error = sqrt((error^2)/newRows);
        end
    end
    answers = newXValues;
    i = iterations;
end

function [x] = gaussian(A,b)
    % GaussPivot: Gauss elimination pivoting
    %   x = GaussPivot(A,b): Gauss elimination with pivoting.
    % input:
    %   A = coefficient matrix
    %   b = right hand side vector
    % output:
    %   x = solution vector
    [m,n]=size(A);
    if m~=n, error('Matrix A must be square'); end
    nb=n+1;
    Aug=[A b];
    % forward elimination
    for k = 1:n-1
      % partial pivoting
      [big,i]=max(abs(Aug(k:n,k)));
      ipr=i+k-1;
      if ipr~=k
        Aug([k,ipr],:)=Aug([ipr,k],:);
      end
      for i = k+1:n
        factor=Aug(i,k)/Aug(k,k);
        Aug(i,k:nb)=Aug(i,k:nb)-factor*Aug(k,k:nb);
      end
    end
    % back substitution
    x=zeros(n,1);
    x(n)=Aug(n,nb)/Aug(n,n);
    for i = n-1:-1:1
      x(i)=(Aug(i,nb)-Aug(i,i+1:n)*x(i+1:n))/Aug(i,i);
    end
end

%used to display the results in the gui
function display_results(roots_results_box)
    global roots;
    roots_results_box.Text = "";
    %roots_results_box.
end


%used to calculate the mean true error
function [mean_t_error] = mean_true_error(A,b,sol,n)
    sum = 0;
    for i = 1:n
        sum = 0;
        for j = 1:n
            sum = sum + (A(i,j) * sol(j));
            row(i) = sum;
        end
    end

    err_row = row;

    err_sum = 0;
    for i = 1:length(row)
        err_arr = abs(err_row(i) - b(i));
        err_sum = err_sum + err_arr;
    end

    mean_t_error = err_sum / length(row);
    


end

%used to specify the stopping criteria (rmse or mae)
function s_criteria(value)
    global stop_c_select;
    stop_c_select = value;
end

%used to specify the method that will be used in calculation
function set_method(m)
    global method;
    method = m.NewValue.Text
end

%used to get the tolerance threshold
function [t] = get_threshold(th_textbox)

    t = str2double(th_textbox.Value);
    disp(t)
    
end

%used to get the starting approximation
function [s_a] = get_s_approx(s_approx_textbox)
    
    arr = s_approx_textbox.Value;
    
    for i = 1 : length(arr)

        %split the entries into string arrays with a comma as the delimeter
        row  = strsplit(string(arr(i)));

        %trim white spaces from all the elements in the row
        row = strtrim(row);
        for j = 1:length(row)
            %convert each string element in the row to a number and assign
            %it a respective place in the matrix
            b(i,j) = str2double(row(j));
            %fprintf("\n %d %d", i, j)
        end
        
    end

    s_a = b;
    
end

%used to add messages to the error box
function [msg] = add_to_message(msg)
    global message;
    message = message + sprintf("\n%s",msg);
end

%used to clear messages from the error box when the calculate button is
%pressed again
function clear_message()
    global message;
    message = "";
end

%check if matrix is diagonally dominant
function [isDiagonalized]= check_matrix(checkMatrix) %Define function
    [rows, cols] = size(checkMatrix); %get the rows and cols
    cols = cols-1; %adjust for not counting col b
    isDiagonallyDom = 1; %assume it is diagonally dominant unless proven otherwise
    notDDRows = []; %a matrix to contain the rows that are not diagonally dominant
    for i = 1: rows %nested for loop to sum each row, and compare the sum to the value of the diagonal number
        rowElements = checkMatrix(i,:);
        sum = 0;
        for j = 1: cols
            if(j~=i)
                sum = sum + abs(rowElements(j));
            end
        end
        if sum >= checkMatrix(i,i)
            %notDDRows(i) = i;
            isDiagonallyDom = 0;
        end
    end
    isDiagonalized = isDiagonallyDom;
    %0 = not diagonolized 1 = diagonalized
    %notDiagonalizedRows = notDDRows;
end


%a function to attempt to make a matrix diagonally dominant
function [new_matrix] = make_dd(aug)
    rows_n = size(aug,1);
    cols_n = size(aug,2);
    
    for count = 1:rows_n %try the algorithm a bunch of times to make sure that rows are moved properly to make diagonal entries the largest entry in each row
        for i = 1:rows_n %cycle through all the rows           
            sum = 0;

            for j = 1:rows_n
                if j~=i
                    sum = sum + abs(aug(i,j));
                end
            end

           if abs(aug(i,i)) < sum 
               temp_row = aug(i,:);
               
                for j = 1:rows_n
                    if j~=i
                        if aug(j,j) > aug(i,i) 
                            aug(i,:) = aug(j,:);
                            aug(j,:) = temp_row;
                        end
                    end
                end

           end 

        end
    end

    new_matrix = aug;
end