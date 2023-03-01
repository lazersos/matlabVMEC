        function ax_keyboard_control(src, event)
            key = event.Key;
            ax = gca;
            %Determine vectors
            position = ax.CameraPosition;
            target = ax.CameraTarget;
            up = ax.CameraUpVector;
            vec = target - position; %Vector from camera position to target
            
            switch key
                case 'w'
                    %move camera forward
                    position = position + 0.1/norm(vec) * vec;
                    target = position + 0.1/norm(vec) * vec;
                    
                case 's'
                    %move camera backward
                    position = position - 0.1/norm(vec) * vec;
                    target = target - 0.1/norm(vec) * vec;
                case 'a'
                    %move camera to the left
                    move = cross(up,vec);
                    position = position + 0.1/norm(move) * move;
                    target = target + 0.1/norm(move) * move;
                    
                case 'd'
                    %move camera to the right
                    move = cross(up,vec);
                    position = position - 0.1/norm(move) * move;
                    target = target - 0.1/norm(move) * move;
                    
                case 'j'
                    %Rotate camera view right (around up axis)
                    M = makehgtform('axisrotate',up,deg2rad(5));
                    vec = M(1:3,1:3) * vec'; %transform sight line vector
                    target = position + vec'; %calculate target "position" vector
                case 'l'
                    M = makehgtform('axisrotate',up,deg2rad(-5)); %negative for inverted rotation
                    vec = M(1:3,1:3) * vec'; %transform sight line vector
                    target = position + vec'; %calculate target "position" vector
                case 'i'
                    %Rotate camera view up
                    M = makehgtform('axisrotate',cross(vec,up),deg2rad(5)); %negative for inverted rotation
                    vec = M(1:3,1:3) * vec'; %transform sight line vector
                    up = M(1:3,1:3) * up';
                    target = position + vec'; %calculate target "position" vector
                case 'k'
                    M = makehgtform('axisrotate',cross(vec,up),deg2rad(-5));
                    vec = M(1:3,1:3) * vec'; %transform sight line vector
                    up = M(1:3,1:3) * up';
                    target = position + vec'; %calculate target "position" vector
                case 'o'
                    M = makehgtform('axisrotate',vec,deg2rad(5));
                    up = M(1:3,1:3) * up';
                case 'u'
                    M = makehgtform('axisrotate',-vec,deg2rad(5));
                    up = M(1:3,1:3) * up';
            end
            ax.CameraUpVector = up;
            ax.CameraPosition = position;
            ax.CameraTarget = target;
        end