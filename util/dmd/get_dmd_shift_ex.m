function the_shift_ex = get_dmd_shift_ex(x_pos, dmd_shift_mode)
% 0<=x_pos<1，实数，对应一圈的不同位置
%  0<=x_pos<1, corresponds to different positions of a circle

x_pos = mod(x_pos, 1);

% 一圈180帧的偏移矫正 Offset correction of 180 frames in a circle
switch dmd_shift_mode 
    case 'none'
        shift_ex = zeros(180, 2);  % 不修正 do not correct
    case 't0'
        shift_ex = -[
-51,-51;-52,-50;-54,-48;-56,-46;-57,-45;-59,-43;-61,-41;-62,-39;-63,-37;-65,-35;-66,-33;-67,-31;-68,-29;-69,-25;-70,-23;-71,-21;-72,-18;-73,-16;-73,-13;-74,-11;-75,-8;-75,-5;-75,-2;-76,0;-76,3;-76,5;-76,7;-76,10;-76,12;-75,15;-75,17;-74,20;-74,22;-73,26;-72,28;-72,30;-71,33;-70,35;-69,37;-68,39;-66,42;-65,44;-64,46;-62,48;-60,50;-59,53;-57,55;-54,57;-52,59;-50,61;-47,63;-45,64;-42,66;-40,68;-37,69;-34,70;-32,72;-29,73;-25,74;-22,75;-19,76;-16,76;-13,77;-10,77;-7,77;-4,78;-1,78;2,78;5,77;8,77;11,77;13,76;16,76;19,75;22,74;25,73;27,72;30,71;33,70;36,68;38,67;40,65;43,64;45,62;47,60;49,59;51,57;53,55;55,53;56,51;58,49;60,47;61,44;62,42;64,40;65,37;66,35;67,32;68,30;69,28;70,25;70,23;71,20;71,17;72,15;72,12;72,10;72,8;72,6;72,3;72,1;72,-1;72,-4;71,-6;71,-8;70,-10;70,-13;69,-15;69,-16;68,-18;67,-20;66,-22;66,-24;65,-27;64,-29;63,-30;62,-32;61,-34;60,-35;59,-37;58,-39;56,-41;55,-42;54,-44;52,-45;51,-47;49,-49;48,-50;46,-52;44,-53;42,-55;40,-56;38,-58;36,-59;34,-60;32,-61;29,-63;26,-64;23,-65;20,-66;17,-67;15,-67;13,-68;10,-68;7,-69;5,-69;2,-69;0,-69;-2,-69;-5,-69;-7,-69;-10,-69;-12,-69;-14,-68;-17,-68;-19,-67;-22,-66;-24,-66;-27,-65;-30,-64;-32,-63;-34,-62;-36,-61;-38,-60;-40,-59;-42,-57;-45,-56;-46,-55;-49,-53;-50,-52;
            ];
    case 't1'
        shift_ex = -[
-52,-52;-53,-51;-55,-49;-57,-47;-58,-46;-60,-44;-62,-42;-63,-39;-64,-37;-67,-35;-68,-33;-69,-31;-70,-29;-71,-26;-72,-24;-73,-22;-74,-19;-75,-17;-75,-13;-76,-11;-77,-8;-77,-5;-77,-2;-78,0;-78,3;-78,5;-78,7;-78,10;-78,13;-77,16;-77,18;-76,21;-76,23;-75,26;-74,28;-74,30;-73,33;-71,35;-70,37;-69,40;-67,43;-66,45;-65,47;-63,49;-61,51;-60,54;-58,56;-55,58;-53,60;-51,62;-48,64;-46,65;-42,67;-40,69;-37,71;-34,72;-32,74;-29,75;-26,76;-23,77;-20,78;-17,78;-14,79;-10,79;-7,79;-4,80;-1,80;2,80;5,79;8,79;11,79;13,78;16,78;20,77;23,76;26,75;28,74;31,73;33,72;36,70;38,69;40,67;43,65;45,63;48,61;50,60;52,58;54,56;56,54;57,52;59,50;61,48;62,45;63,43;65,41;66,38;67,36;68,32;69,30;71,28;72,25;72,23;73,20;73,18;74,16;74,13;74,11;74,9;74,6;74,3;74,1;74,-1;74,-4;73,-6;73,-8;72,-10;72,-13;71,-15;71,-17;70,-19;69,-21;68,-23;68,-25;67,-27;66,-29;65,-30;64,-32;62,-34;61,-35;60,-38;59,-40;57,-42;56,-43;55,-45;53,-46;52,-48;50,-50;49,-51;47,-53;45,-54;43,-56;41,-57;38,-59;36,-60;34,-62;32,-63;29,-65;26,-66;24,-67;21,-68;18,-69;16,-69;14,-70;10,-70;7,-71;5,-71;2,-71;0,-71;-2,-71;-5,-71;-7,-71;-10,-71;-12,-71;-15,-70;-18,-70;-20,-69;-23,-68;-25,-68;-27,-67;-30,-66;-32,-65;-34,-64;-36,-62;-38,-61;-41,-60;-43,-58;-46,-57;-47,-56;-50,-54;-51,-53;
            ];
    case 't1u5'  
        shift_ex = -[
-54,-54;-56,-52;-57,-50;-59,-48;-61,-46;-63,-43;-64,-42;-66,-39;-68,-37;-69,-35;-70,-33;-72,-30;-73,-28;-74,-26;-75,-23;-76,-21;-77,-18;-77,-16;-78,-13;-79,-11;-79,-8;-80,-5;-80,-2;-80,0;-80,3;-80,5;-80,8;-80,11;-80,13;-79,16;-79,18;-78,21;-78,23;-77,26;-76,30;-76,30;-75,33;-74,36;-73,38;-71,41;-69,45;-68,47;-67,48;-65,51;-63,54;-60,56;-58,58;-55,61;-53,64;-50,65;-48,67;-45,69;-42,70;-39,72;-36,73;-33,75;-31,76;-28,77;-25,78;-22,78;-19,79;-16,79;-13,80;-10,81;-7,81;-4,81;-1,81;2,80;5,80;8,80;11,80;13,79;16,78;19,78;21,77;24,76;26,75;29,74;31,73;34,71;36,70;38,68;41,67;43,65;45,64;47,62;49,60;51,58;53,56;55,54;56,52;58,50;60,48;61,46;62,44;64,41;65,39;66,36;67,34;68,31;69,29;70,27;70,24;71,22;72,19;72,16;72,14;73,11;73,8;73,6;73,3;73,1;73,-2;73,-5;72,-7;72,-9;71,-12;71,-15;70,-17;69,-19;69,-20;68,-22;67,-25;66,-27;65,-29;64,-31;62,-33;61,-36;61,-36;59,-39;57,-41;55,-43;55,-43;54,-45;53,-47;51,-47;50,-50;49,-50;47,-52;45,-54;43,-55;41,-57;38,-58;36,-60;34,-62;32,-63;29,-64;26,-66;23,-67;20,-68;17,-69;15,-69;14,-70;10,-71;8,-71;5,-71;2,-71;0,-72;-2,-72;-5,-72;-7,-72;-10,-72;-12,-71;-15,-71;-17,-71;-20,-70;-22,-70;-25,-69;-26,-69;-29,-68;-32,-66;-34,-66;-36,-65;-39,-64;-41,-62;-43,-61;-46,-60;-47,-58;-50,-57;-52,-55;        
            ];
    case 't1u6' 
        shift_ex = -[
-50,-50;-52,-48;-53,-47;-55,-45;-57,-43;-59,-40;-60,-39;-62,-36;-63,-34;-64,-33;-65,-31;-67,-28;-68,-26;-69,-24;-70,-21;-71,-20;-71,-17;-72,-15;-73,-12;-73,-10;-74,-8;-74,-5;-75,-2;-75,0;-75,3;-75,5;-75,7;-75,10;-75,12;-75,15;-74,17;-74,20;-74,21;-73,24;-71,28;-72,28;-71,31;-70,33;-69,36;-67,38;-65,42;-64,44;-63,45;-61,48;-59,50;-57,53;-55,55;-52,57;-50,60;-47,61;-45,63;-42,65;-40,66;-37,68;-34,69;-31,71;-29,71;-26,72;-23,73;-21,74;-18,75;-15,75;-12,75;-9,76;-6,76;-4,76;-1,76;2,76;5,75;8,75;11,75;12,74;15,73;18,73;20,72;23,71;25,70;27,69;29,68;32,66;34,66;36,64;39,62;41,61;42,60;44,58;46,56;48,54;50,52;52,50;53,49;55,47;56,44;57,43;58,41;60,38;61,36;62,33;63,32;64,29;65,27;65,25;66,22;67,20;67,17;68,15;68,13;68,10;68,7;68,5;68,2;68,1;68,-2;67,-5;67,-7;67,-9;66,-12;65,-14;65,-16;64,-18;64,-19;63,-21;62,-23;61,-25;60,-27;59,-29;58,-31;56,-34;56,-33;54,-36;53,-38;52,-41;51,-40;50,-42;49,-43;48,-44;46,-46;46,-46;44,-48;42,-50;40,-51;38,-53;36,-54;34,-56;31,-57;29,-58;27,-60;24,-61;21,-62;18,-63;16,-64;14,-64;13,-65;9,-66;7,-66;5,-66;2,-66;0,-67;-2,-67;-5,-67;-7,-67;-9,-67;-11,-66;-14,-66;-16,-66;-19,-65;-20,-65;-23,-64;-24,-64;-27,-63;-30,-62;-31,-61;-33,-60;-36,-59;-38,-58;-40,-57;-42,-55;-44,-54;-46,-53;-48,-51;
            ];
    case 't1u7'
        shift_ex = -[
-58,-58;-59,-57;-61,-54;-63,-52;-65,-51;-66,-48;-69,-46;-70,-43;-72,-41;-74,-38;-75,-36;-76,-34;-77,-32;-79,-28;-80,-26;-81,-24;-81,-21;-83,-18;-83,-14;-84,-12;-85,-9;-85,-5;-85,-2;-85,0;-85,3;-85,6;-86,8;-85,11;-85,14;-84,18;-84,20;-83,23;-83,25;-82,29;-81,31;-80,33;-79,36;-78,39;-77,41;-76,44;-73,47;-72,49;-71,51;-69,54;-67,56;-65,59;-63,61;-60,64;-58,66;-55,68;-53,70;-50,72;-47,74;-44,76;-41,78;-37,79;-35,80;-32,81;-29,83;-25,84;-22,85;-19,85;-16,86;-11,86;-8,87;-5,87;-2,87;2,87;5,87;8,87;11,86;14,86;17,85;21,84;24,83;28,82;30,81;33,80;35,79;39,77;41,76;43,74;47,72;49,70;52,68;54,66;56,64;59,62;61,60;62,58;64,56;66,53;68,51;69,48;71,46;73,43;74,41;76,37;76,34;77,31;79,28;79,26;80,23;80,21;81,18;81,15;81,13;81,11;82,7;82,4;82,2;82,0;82,-4;81,-6;81,-8;80,-11;80,-14;79,-16;79,-18;78,-21;77,-23;76,-25;75,-27;74,-30;73,-32;72,-33;71,-35;70,-38;69,-40;67,-42;65,-44;64,-47;63,-48;61,-50;59,-52;58,-54;56,-56;54,-57;52,-59;50,-60;48,-63;46,-64;42,-66;40,-68;38,-69;35,-70;32,-72;29,-73;26,-74;23,-75;20,-76;18,-77;15,-77;11,-78;8,-78;6,-79;2,-79;0,-79;-2,-79;-6,-78;-8,-79;-11,-78;-13,-78;-17,-77;-20,-77;-22,-76;-25,-75;-27,-75;-30,-74;-33,-73;-35,-72;-38,-71;-40,-70;-43,-68;-45,-67;-48,-65;-50,-63;-52,-62;-55,-60;-56,-59;
            ];
    case 't1u8'  % 联合矫正
        shift_ex = -[        
-60,-60;-61,-59;-63,-56;-65,-54;-67,-52;-69,-50;-71,-47;-73,-45;-75,-42;-76,-39;-78,-37;-79,-34;-80,-32;-82,-29;-83,-26;-84,-23;-85,-21;-86,-18;-86,-15;-87,-11;-88,-9;-88,-5;-88,-2;-89,0;-89,3;-89,6;-89,9;-88,12;-88,15;-87,19;-87,21;-86,24;-85,26;-85,30;-84,32;-83,34;-82,38;-81,41;-80,43;-78,46;-76,49;-75,51;-74,53;-72,56;-70,59;-68,61;-65,63;-62,67;-60,69;-57,71;-53,74;-51,76;-47,78;-41,80;-41,81;-38,83;-34,85;-32,85;-28,87;-25,87;-21,88;-17,89;-14,89;-11,90;-7,90;-4,91;-1,90;2,90;6,90;9,90;12,90;15,89;19,88;22,87;25,86;29,85;31,84;34,83;36,82;40,80;42,79;45,77;48,75;50,73;53,71;55,69;57,67;59,65;62,63;63,61;65,58;67,56;69,54;71,51;72,49;74,46;75,44;76,41;77,38;79,36;80,32;80,30;81,27;82,24;83,21;83,18;84,15;84,12;84,10;84,6;84,3;84,0;84,-2;84,-5;83,-8;83,-11;82,-14;81,-16;81,-18;80,-21;79,-24;78,-26;78,-28;76,-31;75,-33;75,-34;73,-36;72,-39;71,-41;69,-43;68,-46;66,-48;65,-49;63,-52;61,-54;59,-55;57,-57;56,-59;53,-61;52,-62;49,-64;47,-66;44,-68;41,-70;39,-71;36,-73;33,-74;30,-76;27,-77;24,-78;21,-79;19,-79;15,-79;11,-80;9,-81;6,-81;3,-81;0,-81;-2,-82;-5,-81;-8,-82;-11,-81;-13,-81;-15,-80;-19,-80;-21,-79;-24,-79;-27,-78;-30,-77;-32,-76;-35,-75;-38,-74;-40,-73;-43,-71;-45,-70;-48,-68;-51,-66;-53,-65;-55,-63;-58,-61;
        ];
    case 't1u9'  % 联合矫正
        shift_ex = -[
-57,-57;-58,-56;-62,-52;-63,-51;-65,-50;-66,-47;-68,-45;-69,-42;-71,-40;-73,-38;-74,-35;-75,-34;-77,-31;-78,-27;-79,-26;-80,-24;-80,-20;-81,-17;-82,-13;-83,-12;-83,-8;-84,-5;-84,-2;-83,0;-83,3;-83,6;-84,8;-84,11;-84,14;-83,18;-82,20;-82,23;-82,25;-80,28;-80,31;-79,33;-77,35;-76,38;-75,40;-75,43;-72,46;-71,48;-70,50;-68,53;-65,55;-64,58;-62,60;-59,63;-57,64;-54,67;-52,69;-49,70;-46,73;-44,75;-40,77;-36,77;-34,78;-31,79;-28,81;-24,82;-22,83;-19,83;-16,85;-11,84;-8,85;-5,85;-2,85;2,85;5,85;8,85;11,84;14,84;17,83;21,82;24,81;26,80;29,79;32,78;34,77;37,75;40,74;42,72;45,71;48,69;50,67;52,65;55,63;57,61;59,59;60,57;63,56;64,52;66,50;67,47;69,45;71,42;72,40;74,36;75,33;75,30;77,28;78,25;78,22;78,20;78,18;79,15;79,13;79,11;80,7;80,4;80,2;80,0;80,-4;79,-6;80,-8;78,-11;78,-14;78,-16;77,-18;76,-21;75,-23;75,-25;73,-26;73,-30;72,-31;70,-32;70,-34;69,-37;67,-39;66,-41;63,-43;63,-46;62,-47;60,-49;58,-51;57,-53;55,-55;53,-56;51,-58;50,-59;48,-62;45,-63;41,-65;39,-67;38,-68;34,-68;32,-71;28,-71;25,-72;23,-73;20,-74;18,-75;15,-77;11,-77;8,-77;6,-78;2,-78;0,-79;-2,-78;-6,-77;-8,-77;-11,-77;-13,-77;-17,-76;-19,-76;-22,-76;-25,-74;-26,-74;-29,-73;-33,-72;-34,-71;-37,-70;-39,-69;-42,-68;-44,-66;-47,-65;-49,-63;-51,-61;-54,-60;-55,-59;
        ];
    case 't1u10' % 正面矫正 
        shift_ex = -[
-60,-60;-61,-59;-65,-54;-66,-54;-68,-52;-70,-50;-71,-47;-72,-44;-74,-42;-76,-40;-77,-36;-78,-35;-80,-32;-81,-28;-83,-27;-84,-25;-84,-21;-84,-18;-86,-14;-86,-12;-87,-8;-89,-5;-89,-2;-88,0;-87,3;-86,6;-87,8;-88,12;-88,15;-86,19;-85,21;-86,24;-85,26;-83,29;-84,33;-83,35;-81,37;-79,39;-79,42;-78,45;-75,48;-74,50;-74,53;-71,55;-67,57;-66,60;-65,63;-61,66;-59,67;-56,70;-54,72;-50,72;-48,75;-46,79;-42,81;-38,81;-36,82;-32,83;-29,84;-25,86;-23,86;-20,87;-17,89;-12,88;-9,89;-5,89;-2,89;2,89;5,89;8,89;11,88;14,88;18,87;22,86;25,85;27,84;30,83;33,82;35,81;39,79;41,77;43,75;47,74;50,72;52,70;54,68;57,66;59,64;60,61;62,60;65,58;67,55;69,52;70,49;72,47;74,44;75,42;78,38;78,35;78,32;80,30;82,27;82,24;82,21;81,19;83,16;82,14;82,12;84,8;83,4;83,2;84,0;82,-4;83,-6;83,-8;81,-11;81,-14;81,-16;80,-19;79,-22;79,-24;78,-26;76,-27;77,-31;75,-32;72,-33;73,-35;71,-38;71,-41;68,-42;65,-44;65,-48;64,-49;62,-51;60,-53;59,-55;57,-57;55,-58;53,-60;52,-61;50,-65;47,-65;42,-67;41,-71;39,-70;36,-72;34,-75;30,-75;26,-76;24,-76;21,-76;19,-78;16,-81;11,-80;8,-80;6,-81;2,-82;0,-84;-2,-82;-6,-81;-8,-81;-12,-81;-14,-81;-18,-79;-20,-80;-23,-80;-26,-78;-27,-77;-31,-77;-34,-75;-36,-75;-38,-73;-41,-73;-44,-71;-46,-69;-49,-68;-51,-66;-50,-63;-56,-63;-57,-62;
        ];
    case 't1u11' 
        shift_ex = -[
-65,-65;-66,-64;-70,-59;-72,-58;-73,-56;-75,-54;-77,-51;-79,-48;-81,-46;-83,-43;-84,-40;-86,-39;-87,-35;-89,-31;-89,-29;-90,-27;-92,-23;-93,-19;-93,-15;-94,-13;-95,-9;-95,-6;-95,-3;-95,0;-95,3;-95,7;-95,9;-95,13;-95,16;-94,20;-94,23;-93,26;-92,28;-91,32;-90,35;-90,38;-88,41;-87,44;-86,46;-85,49;-82,53;-80,55;-79,57;-77,61;-74,64;-72,66;-70,68;-67,72;-65,74;-61,76;-59,78;-56,80;-52,83;-49,84;-45,87;-41,88;-39,89;-36,91;-32,92;-28,93;-25,94;-22,95;-18,95;-13,96;-10,96;-6,96;-3,97;1,97;5,96;8,96;12,96;15,95;19,95;23,93;26,92;29,91;32,90;35,89;38,88;41,86;44,84;47,83;50,81;53,79;56,77;58,75;61,72;64,70;66,68;68,66;69,63;72,60;74,58;76,55;78,52;80,49;81,46;83,42;84,39;86,36;87,33;88,29;89,27;89,24;90,22;90,18;91,16;91,14;91,9;91,6;91,3;91,1;91,-4;91,-6;90,-8;90,-12;89,-15;89,-17;88,-20;87,-23;86,-26;85,-28;85,-30;83,-33;82,-35;81,-37;80,-39;78,-42;77,-44;75,-47;73,-49;72,-52;70,-53;68,-56;66,-58;65,-60;62,-62;61,-64;58,-66;57,-67;54,-70;51,-72;47,-75;44,-76;43,-77;39,-79;36,-80;32,-82;29,-83;27,-84;23,-85;21,-86;17,-87;13,-88;10,-88;7,-88;3,-89;0,-89;-2,-89;-7,-89;-9,-89;-12,-88;-15,-88;-19,-87;-22,-87;-25,-86;-28,-85;-30,-85;-33,-84;-37,-82;-39,-81;-42,-80;-45,-79;-48,-77;-50,-76;-53,-74;-56,-72;-58,-70;-61,-68;-63,-67;
        ];
    case 't1v1'  % 侧面好用
        shift_ex = -[
-54,-54;-56,-52;-57,-50;-59,-48;-61,-46;-62,-44;-64,-42;-66,-39;-67,-37;-69,-35;-70,-33;-71,-31;-73,-28;-74,-26;-74,-23;-76,-21;-76,-18;-77,-16;-78,-13;-79,-11;-79,-8;-79,-5;-79,-3;-80,0;-80,3;-80,5;-80,8;-79,11;-80,13;-79,16;-79,18;-78,21;-78,23;-77,26;-75,29;-76,30;-75,33;-73,35;-72,38;-70,40;-68,44;-66,45;-66,47;-63,49;-61,52;-58,54;-56,56;-53,59;-50,61;-48,63;-46,65;-43,66;-40,68;-37,70;-35,72;-32,73;-29,74;-26,75;-24,76;-21,77;-18,78;-15,78;-12,79;-9,80;-6,80;-4,81;-1,81;2,80;5,80;8,80;11,80;13,79;16,79;19,78;21,77;24,76;26,75;29,74;31,73;33,72;36,70;38,68;41,67;42,66;45,64;47,62;49,61;50,59;52,57;54,55;56,52;57,51;59,49;60,47;62,44;63,42;64,40;65,37;66,35;67,32;68,30;69,28;69,25;70,22;71,20;71,17;71,15;71,12;71,9;72,7;71,4;72,2;71,-1;71,-3;70,-6;71,-8;70,-11;69,-13;68,-16;68,-18;68,-19;68,-22;66,-24;66,-27;64,-29;63,-31;62,-33;60,-35;61,-36;58,-39;57,-41;55,-43;55,-43;54,-45;52,-46;51,-47;49,-49;48,-49;46,-51;44,-53;42,-54;40,-56;38,-57;35,-59;33,-60;31,-61;28,-63;25,-65;22,-66;20,-67;17,-69;15,-69;13,-68;9,-70;8,-70;5,-71;2,-71;0,-72;-2,-71;-5,-72;-7,-71;-10,-71;-12,-71;-15,-71;-17,-71;-20,-70;-22,-71;-25,-69;-26,-68;-29,-68;-32,-66;-34,-66;-36,-65;-39,-64;-41,-62;-43,-61;-46,-60;-47,-58;-50,-57;-52,-55;
            ];
    case 't1w1'  % 正面好用
        shift_ex = -[
-52,-52;-54,-50;-56,-49;-58,-47;-60,-45;-62,-42;-63,-41;-65,-38;-66,-36;-67,-34;-69,-32;-70,-29;-71,-27;-72,-25;-73,-22;-74,-20;-75,-18;-76,-16;-76,-13;-77,-11;-78,-8;-78,-5;-78,-2;-78,0;-79,3;-79,5;-79,8;-78,11;-78,13;-78,16;-78,18;-77,21;-77,23;-76,26;-74,29;-75,29;-74,32;-73,35;-72,37;-70,40;-68,44;-67,46;-66,47;-64,50;-62,53;-59,55;-57,57;-54,60;-52,62;-49,64;-47,66;-44,68;-41,69;-38,71;-36,72;-32,73;-30,74;-27,75;-24,76;-22,77;-19,78;-16,78;-13,78;-10,79;-7,79;-4,79;-1,79;2,79;5,79;8,78;11,78;13,77;16,77;19,76;21,75;23,74;25,73;28,72;30,71;33,70;35,68;37,67;40,65;42,64;44,62;46,60;48,59;50,57;52,55;54,53;55,51;57,49;58,47;59,45;61,43;63,40;64,38;65,35;66,33;67,30;67,28;68,26;69,24;70,21;70,18;71,16;71,14;71,11;72,8;72,6;71,3;71,1;71,-2;71,-5;70,-7;70,-9;69,-12;69,-15;68,-17;67,-19;67,-20;66,-22;65,-24;64,-26;63,-28;62,-30;61,-33;59,-35;59,-35;57,-38;56,-40;54,-42;54,-42;52,-44;51,-45;50,-46;48,-48;48,-49;46,-51;44,-53;42,-54;40,-56;38,-58;36,-59;33,-61;31,-62;29,-63;25,-65;22,-66;19,-67;16,-67;15,-68;14,-68;10,-69;8,-70;5,-70;2,-70;0,-70;-2,-70;-5,-70;-7,-70;-10,-70;-12,-69;-15,-69;-17,-69;-20,-68;-21,-68;-24,-67;-25,-67;-28,-66;-31,-65;-33,-64;-35,-63;-38,-62;-40,-61;-42,-60;-45,-58;-46,-57;-49,-55;-51,-53;
            ]; 
    case 't2'
        shift_ex = -[
-53,-53;-54,-52;-56,-50;-58,-48;-59,-47;-61,-45;-63,-43;-64,-40;-65,-38;-68,-36;-69,-34;-70,-32;-71,-30;-72,-26;-73,-24;-74,-22;-75,-19;-76,-17;-76,-13;-77,-11;-78,-8;-78,-5;-78,-2;-79,0;-79,3;-79,5;-79,7;-79,10;-79,13;-78,16;-78,18;-77,21;-77,23;-76,27;-75,29;-75,31;-74,34;-72,36;-71,38;-70,41;-68,44;-67,46;-66,48;-64,50;-62,52;-61,55;-59,57;-56,59;-54,61;-52,63;-49,65;-47,66;-43,68;-41,70;-38,72;-35,73;-33,75;-30,76;-26,77;-23,78;-20,79;-17,79;-14,80;-10,80;-7,80;-4,81;-1,81;2,81;5,80;8,80;11,80;13,79;16,79;20,78;23,77;26,76;28,75;31,74;34,73;37,71;39,70;41,68;44,66;46,64;49,62;51,61;53,59;55,57;57,55;58,53;60,51;62,49;63,46;64,44;66,42;67,39;68,37;69,33;70,31;72,29;73,26;73,24;74,21;74,18;75,16;75,13;75,11;75,9;75,6;75,3;75,1;75,-1;75,-4;74,-6;74,-8;73,-10;73,-13;72,-15;72,-17;71,-19;70,-21;69,-23;69,-25;68,-28;67,-30;66,-31;65,-33;63,-35;62,-36;61,-39;60,-41;58,-43;57,-44;56,-46;54,-47;53,-49;51,-51;50,-52;48,-54;46,-55;44,-57;42,-58;39,-60;37,-61;35,-63;33,-64;30,-66;27,-67;24,-68;21,-69;18,-70;16,-70;14,-71;10,-71;7,-72;5,-72;2,-72;0,-72;-2,-72;-5,-72;-7,-72;-10,-72;-12,-72;-15,-71;-18,-71;-20,-70;-23,-69;-25,-69;-28,-68;-31,-67;-33,-66;-35,-65;-37,-63;-39,-62;-42,-61;-44,-59;-47,-58;-48,-57;-51,-55;-52,-54;
            ];
    case 't3' % 正面
        shift_ex = -[
-65,-66;-67,-64;-70,-60;-72,-59;-74,-55;-75,-54;-78,-51;-79,-47;-81,-46;-82,-43;-84,-40;-85,-38;-86,-33;-87,-31;-89,-28;-90,-24;-90,-22;-91,-18;-92,-15;-92,-13;-93,-9;-93,-5;-93,-4;-93,0;-93,4;-93,6;-93,10;-92,14;-92,17;-91,20;-90,23;-89,26;-88,31;-87,32;-86,35;-85,38;-84,41;-82,44;-81,47;-79,49;-77,52;-75,55;-73,58;-72,60;-69,62;-67,65;-65,68;-63,69;-61,71;-58,74;-56,75;-52,77;-50,79;-48,80;-44,82;-41,83;-39,84;-36,85;-33,87;-30,87;-27,89;-25,89;-21,90;-18,91;-15,91;-11,92;-9,92;-5,92;-2,92;1,92;5,92;9,92;12,91;15,91;17,90;21,90;23,89;27,88;30,87;32,86;36,84;39,83;41,81;45,80;46,79;50,77;53,75;55,73;57,71;59,69;62,67;64,65;66,63;69,60;70,58;73,55;75,52;75,50;78,47;79,45;81,41;82,40;84,36;84,34;86,31;87,28;88,24;88,22;89,18;90,14;90,12;91,8;91,6;91,1;91,-2;91,-5;91,-8;90,-11;90,-14;89,-18;89,-21;88,-26;87,-27;85,-32;85,-34;83,-36;83,-38;81,-42;79,-44;78,-47;76,-50;74,-53;73,-55;71,-58;68,-60;67,-62;65,-65;62,-67;60,-69;58,-71;55,-73;52,-75;49,-77;45,-79;44,-80;40,-82;38,-83;34,-85;31,-86;28,-87;26,-88;22,-89;18,-90;15,-90;12,-91;9,-91;7,-92;1,-91;0,-92;-4,-92;-7,-91;-11,-91;-14,-91;-17,-90;-20,-90;-23,-89;-26,-88;-30,-87;-32,-86;-36,-85;-39,-83;-42,-82;-45,-80;-47,-79;-50,-78;-53,-76;-56,-74;-58,-72;-61,-70;-63,-68;
            ];
    case 't3u1' % 正侧面联合
        shift_ex = -[
-63,-68;-66,-66;-68,-64;-71,-62;-73,-59;-75,-57;-77,-54;-79,-51;-81,-48;-83,-46;-84,-43;-86,-40;-87,-37;-89,-34;-90,-31;-91,-27;-92,-24;-93,-20;-94,-17;-95,-14;-95,-10;-96,-7;-96,-4;-96,0;-96,4;-96,6;-96,10;-95,14;-95,18;-94,21;-94,24;-93,27;-91,32;-91,33;-89,36;-88,39;-87,42;-85,46;-84,49;-82,51;-80,54;-78,57;-76,60;-74,62;-72,65;-69,67;-67,70;-65,72;-63,74;-60,76;-58,77;-54,80;-52,82;-50,83;-45,85;-43,86;-40,87;-37,88;-34,90;-31,90;-28,92;-26,92;-22,93;-18,94;-15,94;-11,94;-7,95;-4,95;0,95;5,95;7,94;10,94;14,93;18,92;21,91;24,91;28,90;30,88;34,87;37,86;40,84;42,83;45,81;48,79;51,77;53,75;56,74;58,71;61,69;63,67;65,65;67,62;69,60;71,57;73,55;74,52;76,49;77,47;79,44;81,41;82,38;83,36;84,33;85,30;86,27;86,25;87,22;87,19;88,16;88,13;89,10;89,7;89,5;89,1;89,-2;88,-5;89,-8;88,-11;88,-14;87,-18;86,-20;85,-25;84,-26;82,-31;82,-33;81,-35;80,-37;78,-40;77,-43;75,-45;73,-48;71,-51;70,-53;68,-56;66,-58;64,-59;62,-62;60,-64;58,-67;56,-69;53,-70;50,-72;47,-75;44,-77;43,-78;39,-80;37,-81;33,-83;30,-84;27,-85;25,-85;21,-87;18,-88;15,-88;12,-89;9,-89;7,-90;2,-90;0,-90;-4,-90;-7,-90;-10,-90;-13,-90;-16,-90;-19,-89;-22,-88;-25,-88;-27,-87;-30,-86;-33,-85;-36,-84;-39,-83;-42,-82;-45,-80;-47,-79;-50,-78;-53,-76;-56,-74;-58,-72;-61,-71;
        ];
    otherwise
        fprintf('[Info] get_dmd_shift_ex: use [t1u5] shift \n');
        shift_ex = -[
-54,-54;-56,-52;-57,-50;-59,-48;-61,-46;-63,-43;-64,-42;-66,-39;-68,-37;-69,-35;-70,-33;-72,-30;-73,-28;-74,-26;-75,-23;-76,-21;-77,-18;-77,-16;-78,-13;-79,-11;-79,-8;-80,-5;-80,-2;-80,0;-80,3;-80,5;-80,8;-80,11;-80,13;-79,16;-79,18;-78,21;-78,23;-77,26;-76,30;-76,30;-75,33;-74,36;-73,38;-71,41;-69,45;-68,47;-67,48;-65,51;-63,54;-60,56;-58,58;-55,61;-53,64;-50,65;-48,67;-45,69;-42,70;-39,72;-36,73;-33,75;-31,76;-28,77;-25,78;-22,78;-19,79;-16,79;-13,80;-10,81;-7,81;-4,81;-1,81;2,80;5,80;8,80;11,80;13,79;16,78;19,78;21,77;24,76;26,75;29,74;31,73;34,71;36,70;38,68;41,67;43,65;45,64;47,62;49,60;51,58;53,56;55,54;56,52;58,50;60,48;61,46;62,44;64,41;65,39;66,36;67,34;68,31;69,29;70,27;70,24;71,22;72,19;72,16;72,14;73,11;73,8;73,6;73,3;73,1;73,-2;73,-5;72,-7;72,-9;71,-12;71,-15;70,-17;69,-19;69,-20;68,-22;67,-25;66,-27;65,-29;64,-31;62,-33;61,-36;61,-36;59,-39;57,-41;55,-43;55,-43;54,-45;53,-47;51,-47;50,-50;49,-50;47,-52;45,-54;43,-55;41,-57;38,-58;36,-60;34,-62;32,-63;29,-64;26,-66;23,-67;20,-68;17,-69;15,-69;14,-70;10,-71;8,-71;5,-71;2,-71;0,-72;-2,-72;-5,-72;-7,-72;-10,-72;-12,-71;-15,-71;-17,-71;-20,-70;-22,-70;-25,-69;-26,-69;-29,-68;-32,-66;-34,-66;-36,-65;-39,-64;-41,-62;-43,-61;-46,-60;-47,-58;-50,-57;-52,-55;        
            ];
end

% shift_ex = -[
% %-50,-48;-52,-46;-52,-44;-53,-43;-54,-40;-55,-38;-56,-36;-57,-34;-58,-31;-59,-29;-60,-27;-61,-24;-62,-22;-62,-20;-63,-18;-63,-16;-63,-14;-63,-12;-62,-10;-62,-8;-62,-6;-62,-3;-62,-1;-62,1;-61,3;-61,5;-60,8;-59,10;-59,12;-58,14;-57,16;-56,18;-55,20;-54,22;-53,23;-52,25;-51,27;-50,28;-49,29;-48,31;-46,32;-45,34;-44,35;-43,36;-42,37;-41,38;-40,39;-39,40;-37,41;-35,42;-34,43;-33,44;-31,45;-30,45;-29,46;-27,47;-26,47;-24,48;-22,49;-21,50;-20,51;-19,51;-17,52;-16,53;-14,53;-12,54;-10,54;-9,54;-7,55;-5,55;-4,55;-2,55;0,55;2,55;4,55;6,56;7,56;9,56;11,57;12,57;14,56;16,56;17,56;19,55;21,55;24,55;26,54;28,53;30,53;32,51;34,51;36,50;38,49;40,48;42,47;44,45;45,44;47,43;49,42;51,40;52,40;54,38;56,36;57,34;58,32;60,29;61,27;63,25;64,23;65,21;66,18;67,17;68,14;69,11;69,9;70,7;70,4;71,2;71,0;71,-3;72,-6;72,-9;72,-12;71,-15;70,-18;70,-21;69,-23;69,-25;68,-27;68,-30;67,-32;66,-35;64,-38;63,-41;61,-44;59,-47;57,-49;56,-51;54,-53;52,-55;50,-57;48,-59;45,-61;42,-63;40,-65;38,-66;35,-67;32,-68;30,-69;27,-70;25,-71;22,-72;19,-73;16,-73;12,-74;10,-74;7,-75;4,-75;1,-75;-2,-74;-5,-74;-8,-74;-11,-73;-13,-73;-16,-72;-18,-71;-20,-70;-23,-70;-25,-68;-28,-67;-30,-65;-33,-64;-36,-62;-38,-61;-40,-59;-42,-57;-44,-55;-46,-53;-48,-51;-49,-50;
% %-50,-59;-52,-57;-54,-56;-56,-54;-58,-53;-60,-50;-62,-48;-63,-46;-65,-43;-67,-41;-68,-38;-69,-36;-70,-33;-71,-31;-72,-28;-73,-26;-74,-23;-75,-21;-76,-18;-77,-15;-77,-12;-77,-9;-77,-5;-77,-2;-77,1;-77,3;-77,6;-76,9;-76,12;-76,14;-75,17;-74,20;-74,23;-72,26;-72,28;-70,31;-69,33;-67,35;-65,37;-63,40;-62,41;-60,43;-59,46;-57,48;-56,50;-54,52;-52,53;-49,55;-48,56;-45,58;-43,59;-41,60;-39,62;-36,63;-33,65;-31,66;-29,66;-27,67;-24,67;-22,68;-20,68;-17,69;-14,69;-12,70;-10,70;-8,70;-6,70;-4,70;-2,70;0,70;3,70;5,69;8,69;10,68;13,67;15,67;17,66;20,66;22,65;24,64;26,63;28,62;30,61;32,60;34,59;36,58;38,57;40,56;41,55;43,53;44,51;46,49;47,48;49,46;50,44;52,43;53,41;54,40;56,38;57,37;59,34;60,32;61,30;61,27;62,25;63,24;64,22;65,20;66,18;67,16;67,13;67,11;68,9;68,6;68,4;69,2;69,0;69,-2;69,-4;69,-7;68,-9;67,-12;67,-15;67,-17;67,-20;67,-22;66,-24;66,-26;65,-28;64,-30;63,-33;62,-35;61,-38;60,-40;59,-42;57,-44;56,-47;54,-49;52,-51;51,-53;49,-55;48,-56;46,-58;45,-60;42,-62;40,-63;38,-65;35,-66;33,-68;31,-69;28,-70;25,-71;22,-72;20,-72;17,-73;15,-74;12,-74;10,-75;7,-76;4,-76;2,-76;-1,-76;-4,-76;-7,-76;-10,-76;-13,-76;-16,-76;-19,-76;-21,-75;-24,-74;-26,-73;-30,-72;-32,-71;-34,-69;-37,-68;-40,-66;-42,-65;-44,-64;-47,-62;-49,-60;
% %-34,-64;-36,-63;-38,-62;-40,-61;-42,-59;-44,-58;-46,-56;-48,-55;-50,-53;-52,-51;-54,-49;-56,-48;-57,-46;-59,-44;-61,-41;-62,-39;-64,-37;-65,-35;-67,-32;-68,-30;-69,-27;-70,-25;-71,-22;-72,-20;-73,-17;-74,-14;-74,-12;-75,-9;-75,-6;-76,-4;-76,-1;-76,2;-76,4;-76,7;-75,10;-75,13;-75,15;-74,18;-73,20;-73,23;-72,26;-71,28;-70,31;-69,33;-68,36;-66,38;-65,41;-63,43;-62,45;-60,48;-58,50;-56,52;-55,54;-53,56;-51,58;-48,59;-46,61;-44,63;-42,64;-39,66;-37,67;-34,68;-32,69;-30,70;-27,71;-24,72;-22,73;-19,74;-16,74;-14,75;-11,75;-8,75;-6,76;-3,76;0,76;2,76;5,75;8,75;10,75;13,74;16,74;18,73;21,72;23,71;25,70;28,69;30,68;33,67;35,66;37,64;39,63;42,61;44,60;46,58;48,56;50,55;52,53;53,51;55,49;56,47;58,45;59,43;61,41;62,39;63,37;64,35;66,32;67,30;67,28;68,26;69,23;70,21;70,19;71,16;71,14;71,11;72,9;72,7;72,4;72,2;72,-1;72,-3;71,-6;71,-8;71,-11;70,-13;70,-15;69,-18;68,-20;67,-22;67,-25;66,-27;65,-29;64,-32;62,-34;61,-36;60,-38;58,-40;57,-42;56,-44;54,-46;52,-47;51,-49;49,-51;48,-53;46,-54;44,-56;42,-57;40,-59;38,-60;36,-62;33,-63;31,-64;29,-65;26,-66;24,-67;22,-68;20,-68;17,-69;15,-70;12,-70;10,-71;7,-71;5,-71;2,-71;0,-72;-3,-72;-5,-71;-8,-71;-10,-71;-13,-71;-15,-71;-18,-70;-20,-70;-22,-69;-25,-68;-27,-68;-29,-67;-31,-66;-33,-65;
% %-52,-46;-53,-46;-55,-44;-56,-42;-58,-40;-59,-38;-61,-36;-62,-34;-63,-32;-64,-30;-65,-28;-66,-25;-67,-23;-68,-21;-69,-18;-70,-16;-70,-13;-71,-11;-71,-8;-72,-6;-72,-3;-72,-1;-72,2;-72,5;-72,7;-72,10;-72,12;-72,15;-71,17;-71,20;-70,22;-69,25;-68,28;-67,30;-67,32;-65,35;-64,37;-63,39;-62,41;-60,44;-58,46;-57,48;-55,50;-53,53;-51,54;-49,56;-46,58;-45,60;-42,61;-40,63;-38,65;-36,66;-33,67;-31,68;-28,69;-26,70;-23,71;-21,72;-18,73;-15,73;-13,74;-10,74;-8,74;-5,75;-2,75;1,75;4,75;6,74;9,74;12,74;14,73;16,73;19,72;22,71;24,70;27,69;29,68;31,67;34,66;37,64;39,63;41,61;43,60;45,58;47,57;49,55;51,53;53,51;54,49;56,47;57,45;59,43;60,41;62,39;63,37;64,35;65,32;66,30;67,28;68,25;69,23;69,20;70,18;70,15;71,13;71,11;71,8;71,6;71,3;71,1;71,-2;71,-4;71,-6;70,-8;70,-11;69,-13;69,-16;68,-18;67,-21;66,-23;65,-25;64,-27;63,-29;62,-32;61,-34;59,-36;58,-38;57,-40;55,-42;54,-44;52,-45;50,-47;49,-49;47,-50;45,-52;43,-54;41,-55;40,-56;38,-57;36,-58;34,-60;31,-61;29,-62;27,-63;25,-64;23,-64;21,-65;19,-66;16,-66;14,-67;12,-67;9,-68;7,-68;5,-68;3,-68;0,-68;-2,-68;-5,-68;-7,-68;-9,-68;-12,-67;-14,-67;-16,-66;-19,-66;-22,-65;-24,-64;-26,-63;-28,-62;-31,-62;-32,-61;-35,-60;-37,-58;-39,-57;-41,-56;-43,-55;-44,-53;-46,-52;-48,-50;-50,-49;-51,-47;
% %-53,-45;-54,-45;-56,-43;-57,-41;-59,-39;-60,-37;-62,-35;-63,-32;-64,-30;-65,-28;-66,-26;-67,-23;-68,-21;-69,-19;-69,-16;-70,-14;-70,-11;-71,-9;-71,-6;-72,-4;-72,-1;-72,1;-72,4;-72,7;-72,9;-71,12;-71,14;-71,17;-70,19;-70,22;-69,24;-68,27;-67,30;-66,32;-65,34;-63,37;-62,38;-61,40;-60,43;-58,45;-56,47;-55,49;-53,51;-51,54;-49,55;-47,57;-44,59;-43,60;-40,61;-38,63;-36,65;-34,66;-31,67;-29,68;-26,69;-24,70;-21,71;-19,72;-17,72;-14,72;-12,73;-9,73;-7,73;-4,74;-1,74;2,74;4,74;6,73;9,73;12,73;14,72;16,72;19,71;22,70;24,69;27,68;29,67;31,66;34,65;37,63;39,62;41,60;43,59;45,57;47,56;49,54;51,52;53,50;54,48;56,46;57,44;59,42;60,40;62,38;63,36;64,34;65,31;66,29;67,27;68,24;69,22;69,19;70,17;70,14;71,12;71,10;71,7;71,5;71,2;70,0;70,-3;70,-5;70,-7;69,-9;69,-12;68,-14;68,-17;67,-19;66,-22;65,-24;64,-26;63,-28;62,-30;61,-33;60,-35;58,-37;57,-39;56,-41;54,-42;53,-44;51,-45;49,-47;48,-49;46,-50;44,-52;42,-54;40,-55;39,-56;37,-57;35,-58;33,-60;30,-61;28,-62;26,-63;24,-64;22,-64;20,-65;18,-66;15,-66;13,-67;11,-67;8,-68;6,-68;4,-68;2,-68;-1,-68;-3,-68;-6,-68;-8,-68;-10,-68;-13,-67;-15,-67;-17,-66;-20,-66;-23,-65;-25,-64;-27,-63;-29,-62;-32,-62;-33,-61;-36,-60;-38,-58;-40,-57;-42,-56;-44,-55;-45,-53;-47,-52;-49,-50;-51,-49;-52,-46;
% %-56,-46;-57,-46;-59,-44;-60,-42;-62,-40;-63,-37;-65,-35;-66,-32;-67,-30;-68,-27;-69,-25;-70,-22;-71,-20;-72,-18;-72,-15;-73,-13;-73,-10;-74,-8;-74,-5;-75,-3;-75,0;-75,3;-75,6;-74,9;-74,11;-73,14;-73,16;-73,19;-72,21;-72,25;-70,27;-69,30;-68,33;-67,35;-66,37;-64,40;-62,41;-61,43;-60,46;-58,48;-56,50;-55,52;-53,54;-51,57;-49,58;-46,59;-43,61;-42,62;-39,63;-37,65;-35,67;-33,68;-30,69;-28,70;-25,71;-23,72;-20,73;-18,74;-16,74;-13,73;-11,74;-8,74;-6,74;-3,75;0,75;3,75;5,75;7,73;10,73;13,73;15,72;17,72;20,71;23,70;25,69;28,68;30,67;32,66;35,65;38,63;40,62;42,60;44,58;46,56;47,55;49,53;51,51;53,49;54,47;56,45;57,43;59,41;60,39;62,37;63,35;64,33;65,30;66,28;67,26;68,23;69,21;69,18;70,16;70,13;71,11;71,9;71,6;71,4;71,1;70,-1;70,-4;70,-6;70,-8;69,-10;69,-13;68,-15;68,-18;67,-20;66,-23;65,-25;64,-27;63,-29;62,-31;61,-34;60,-36;58,-38;57,-40;56,-42;54,-43;53,-45;51,-46;49,-48;49,-50;47,-51;45,-53;43,-55;41,-56;40,-57;38,-58;36,-59;34,-61;31,-62;29,-63;27,-64;25,-65;23,-65;21,-67;19,-68;16,-68;14,-69;12,-69;9,-70;7,-70;5,-70;3,-70;0,-70;-2,-70;-5,-71;-7,-71;-9,-71;-12,-70;-14,-70;-17,-69;-20,-69;-23,-68;-25,-67;-27,-66;-29,-65;-33,-64;-34,-63;-37,-62;-39,-60;-42,-59;-44,-58;-46,-57;-47,-55;-49,-54;-51,-52;-53,-51;-54,-48;
% %-56,-48;-57,-46;-59,-44;-60,-42;-62,-40;-63,-38;-65,-35;-66,-33;-68,-30;-69,-28;-70,-25;-71,-23;-72,-20;-73,-17;-73,-15;-74,-12;-74,-10;-75,-7;-75,-4;-75,-1;-75,2;-75,4;-75,7;-75,10;-74,13;-74,15;-73,18;-73,21;-72,23;-71,26;-70,28;-69,31;-68,33;-67,36;-65,39;-64,41;-62,43;-61,45;-59,47;-57,49;-55,51;-53,53;-51,55;-49,57;-47,58;-45,60;-43,61;-41,62;-39,64;-37,65;-34,66;-32,67;-30,68;-28,69;-25,70;-23,70;-20,71;-18,72;-15,72;-13,73;-10,73;-8,73;-6,73;-3,73;-1,73;2,73;4,73;7,72;10,72;12,71;14,70;17,70;19,69;22,68;24,67;26,66;29,65;31,64;33,63;35,61;38,60;40,58;42,56;44,55;46,53;47,51;49,50;51,48;52,46;54,44;55,42;56,40;58,38;59,36;60,34;61,32;62,30;63,27;64,25;64,23;65,20;66,18;67,15;67,13;67,11;68,8;68,6;68,4;68,1;68,-1;68,-4;67,-6;67,-8;67,-10;67,-12;66,-15;65,-17;65,-20;64,-22;63,-24;62,-25;61,-27;61,-29;59,-31;58,-34;57,-36;56,-38;55,-40;53,-42;52,-43;51,-45;49,-46;48,-48;46,-50;45,-51;43,-53;42,-54;40,-56;39,-57;37,-58;35,-59;33,-61;31,-62;29,-63;27,-64;25,-65;23,-66;21,-66;19,-67;17,-68;14,-69;12,-70;10,-70;7,-70;5,-70;2,-71;0,-71;-2,-71;-5,-71;-7,-71;-10,-71;-12,-71;-15,-71;-17,-70;-20,-70;-22,-69;-25,-68;-27,-67;-30,-66;-32,-65;-35,-64;-37,-63;-39,-61;-41,-60;-44,-58;-46,-57;-48,-56;-50,-54;-52,-52;-54,-50;
% %-50,-50;-51,-49;-53,-47;-55,-45;-56,-44;-58,-42;-60,-40;-61,-38;-62,-36;-64,-34;-65,-32;-66,-30;-67,-28;-68,-25;-69,-23;-70,-21;-71,-18;-72,-16;-72,-13;-73,-11;-74,-8;-74,-6;-74,-2;-75,0;-75,3;-75,5;-75,7;-75,10;-75,12;-74,15;-74,17;-73,20;-73,22;-72,25;-71,27;-71,29;-70,32;-69,34;-68,36;-67,38;-65,41;-64,43;-63,45;-61,47;-59,49;-58,52;-56,54;-53,56;-51,58;-49,60;-46,62;-44,63;-41,65;-39,67;-36,68;-33,69;-31,71;-28,72;-25,73;-22,74;-19,75;-16,75;-14,76;-10,76;-7,76;-4,77;-1,77;2,77;5,76;8,76;11,76;13,75;16,75;19,74;22,73;25,72;27,71;30,70;32,69;35,67;37,66;39,64;42,63;44,61;46,59;48,58;50,56;52,54;54,52;55,50;57,48;59,46;60,43;61,41;63,39;64,36;65,34;66,31;67,29;68,27;69,24;69,22;70,19;70,17;71,15;71,12;71,10;71,8;71,6;71,3;71,1;71,-1;71,-4;70,-6;70,-8;69,-10;69,-13;68,-15;68,-16;67,-18;66,-20;65,-22;65,-24;64,-26;63,-28;62,-29;61,-31;60,-33;59,-34;58,-36;57,-38;55,-40;54,-41;53,-43;51,-44;50,-46;48,-48;47,-49;45,-51;43,-52;41,-54;39,-55;37,-57;35,-58;33,-59;31,-60;28,-62;25,-63;23,-64;20,-65;17,-66;15,-66;13,-67;10,-67;7,-68;5,-68;2,-68;0,-68;-2,-68;-5,-68;-7,-68;-10,-68;-12,-68;-14,-67;-17,-67;-19,-66;-22,-65;-24,-65;-26,-64;-29,-63;-31,-62;-33,-61;-35,-60;-37,-59;-39,-58;-41,-56;-44,-55;-45,-54;-48,-52;-49,-51;
% %-50,-50;-51,-49;-53,-47;-55,-45;-56,-44;-58,-42;-60,-40;-61,-38;-62,-36;-64,-34;-65,-32;-66,-30;-67,-28;-68,-25;-69,-23;-70,-21;-71,-18;-72,-16;-72,-13;-73,-11;-74,-8;-74,-5;-74,-2;-75,0;-75,3;-75,5;-75,7;-75,10;-75,12;-74,15;-74,17;-73,20;-73,22;-72,25;-71,27;-71,29;-70,32;-69,34;-68,36;-67,38;-65,41;-64,43;-63,45;-61,47;-59,49;-58,52;-56,54;-53,56;-51,58;-49,60;-46,62;-44,63;-41,65;-39,67;-36,68;-33,69;-31,71;-28,72;-25,73;-22,74;-19,75;-16,75;-13,76;-10,76;-7,76;-4,77;-1,77;2,77;5,76;8,76;11,76;13,75;16,75;19,74;22,73;25,72;27,71;30,70;32,69;35,67;37,66;39,64;42,63;44,61;46,59;48,58;50,56;52,54;54,52;55,50;57,48;59,46;60,43;61,41;63,39;64,36;65,34;66,31;67,29;68,27;69,24;69,22;70,19;70,17;71,15;71,12;71,10;71,8;71,6;71,3;71,1;71,-1;71,-4;70,-6;70,-8;69,-10;69,-13;68,-15;68,-16;67,-18;66,-20;65,-22;65,-24;64,-26;63,-28;62,-29;61,-31;60,-33;59,-34;58,-36;57,-38;55,-40;54,-41;53,-43;51,-44;50,-46;48,-48;47,-49;45,-51;43,-52;41,-54;39,-55;37,-57;35,-58;33,-59;31,-60;28,-62;25,-63;23,-64;20,-65;17,-66;15,-66;13,-67;10,-67;7,-68;5,-68;2,-68;0,-68;-2,-68;-5,-68;-7,-68;-10,-68;-12,-68;-14,-67;-17,-67;-19,-66;-22,-65;-24,-65;-26,-64;-29,-63;-31,-62;-33,-61;-35,-60;-37,-59;-39,-58;-41,-56;-44,-55;-45,-54;-48,-52;-49,-51;
% 
% %t0,t1,t2
% %-51,-51;-52,-50;-54,-48;-56,-46;-57,-45;-59,-43;-61,-41;-62,-39;-63,-37;-65,-35;-66,-33;-67,-31;-68,-29;-69,-25;-70,-23;-71,-21;-72,-18;-73,-16;-73,-13;-74,-11;-75,-8;-75,-5;-75,-2;-76,0;-76,3;-76,5;-76,7;-76,10;-76,12;-75,15;-75,17;-74,20;-74,22;-73,26;-72,28;-72,30;-71,33;-70,35;-69,37;-68,39;-66,42;-65,44;-64,46;-62,48;-60,50;-59,53;-57,55;-54,57;-52,59;-50,61;-47,63;-45,64;-42,66;-40,68;-37,69;-34,70;-32,72;-29,73;-25,74;-22,75;-19,76;-16,76;-13,77;-10,77;-7,77;-4,78;-1,78;2,78;5,77;8,77;11,77;13,76;16,76;19,75;22,74;25,73;27,72;30,71;33,70;36,68;38,67;40,65;43,64;45,62;47,60;49,59;51,57;53,55;55,53;56,51;58,49;60,47;61,44;62,42;64,40;65,37;66,35;67,32;68,30;69,28;70,25;70,23;71,20;71,17;72,15;72,12;72,10;72,8;72,6;72,3;72,1;72,-1;72,-4;71,-6;71,-8;70,-10;70,-13;69,-15;69,-16;68,-18;67,-20;66,-22;66,-24;65,-27;64,-29;63,-30;62,-32;61,-34;60,-35;59,-37;58,-39;56,-41;55,-42;54,-44;52,-45;51,-47;49,-49;48,-50;46,-52;44,-53;42,-55;40,-56;38,-58;36,-59;34,-60;32,-61;29,-63;26,-64;23,-65;20,-66;17,-67;15,-67;13,-68;10,-68;7,-69;5,-69;2,-69;0,-69;-2,-69;-5,-69;-7,-69;-10,-69;-12,-69;-14,-68;-17,-68;-19,-67;-22,-66;-24,-66;-27,-65;-30,-64;-32,-63;-34,-62;-36,-61;-38,-60;-40,-59;-42,-57;-45,-56;-46,-55;-49,-53;-50,-52;
% %-52,-52;-53,-51;-55,-49;-57,-47;-58,-46;-60,-44;-62,-42;-63,-39;-64,-37;-67,-35;-68,-33;-69,-31;-70,-29;-71,-26;-72,-24;-73,-22;-74,-19;-75,-17;-75,-13;-76,-11;-77,-8;-77,-5;-77,-2;-78,0;-78,3;-78,5;-78,7;-78,10;-78,13;-77,16;-77,18;-76,21;-76,23;-75,26;-74,28;-74,30;-73,33;-71,35;-70,37;-69,40;-67,43;-66,45;-65,47;-63,49;-61,51;-60,54;-58,56;-55,58;-53,60;-51,62;-48,64;-46,65;-42,67;-40,69;-37,71;-34,72;-32,74;-29,75;-26,76;-23,77;-20,78;-17,78;-14,79;-10,79;-7,79;-4,80;-1,80;2,80;5,79;8,79;11,79;13,78;16,78;20,77;23,76;26,75;28,74;31,73;33,72;36,70;38,69;40,67;43,65;45,63;48,61;50,60;52,58;54,56;56,54;57,52;59,50;61,48;62,45;63,43;65,41;66,38;67,36;68,32;69,30;71,28;72,25;72,23;73,20;73,18;74,16;74,13;74,11;74,9;74,6;74,3;74,1;74,-1;74,-4;73,-6;73,-8;72,-10;72,-13;71,-15;71,-17;70,-19;69,-21;68,-23;68,-25;67,-27;66,-29;65,-30;64,-32;62,-34;61,-35;60,-38;59,-40;57,-42;56,-43;55,-45;53,-46;52,-48;50,-50;49,-51;47,-53;45,-54;43,-56;41,-57;38,-59;36,-60;34,-62;32,-63;29,-65;26,-66;24,-67;21,-68;18,-69;16,-69;14,-70;10,-70;7,-71;5,-71;2,-71;0,-71;-2,-71;-5,-71;-7,-71;-10,-71;-12,-71;-15,-70;-18,-70;-20,-69;-23,-68;-25,-68;-27,-67;-30,-66;-32,-65;-34,-64;-36,-62;-38,-61;-41,-60;-43,-58;-46,-57;-47,-56;-50,-54;-51,-53;
% -53,-53;-54,-52;-56,-50;-58,-48;-59,-47;-61,-45;-63,-43;-64,-40;-65,-38;-68,-36;-69,-34;-70,-32;-71,-30;-72,-26;-73,-24;-74,-22;-75,-19;-76,-17;-76,-13;-77,-11;-78,-8;-78,-5;-78,-2;-79,0;-79,3;-79,5;-79,7;-79,10;-79,13;-78,16;-78,18;-77,21;-77,23;-76,27;-75,29;-75,31;-74,34;-72,36;-71,38;-70,41;-68,44;-67,46;-66,48;-64,50;-62,52;-61,55;-59,57;-56,59;-54,61;-52,63;-49,65;-47,66;-43,68;-41,70;-38,72;-35,73;-33,75;-30,76;-26,77;-23,78;-20,79;-17,79;-14,80;-10,80;-7,80;-4,81;-1,81;2,81;5,80;8,80;11,80;13,79;16,79;20,78;23,77;26,76;28,75;31,74;34,73;37,71;39,70;41,68;44,66;46,64;49,62;51,61;53,59;55,57;57,55;58,53;60,51;62,49;63,46;64,44;66,42;67,39;68,37;69,33;70,31;72,29;73,26;73,24;74,21;74,18;75,16;75,13;75,11;75,9;75,6;75,3;75,1;75,-1;75,-4;74,-6;74,-8;73,-10;73,-13;72,-15;72,-17;71,-19;70,-21;69,-23;69,-25;68,-28;67,-30;66,-31;65,-33;63,-35;62,-36;61,-39;60,-41;58,-43;57,-44;56,-46;54,-47;53,-49;51,-51;50,-52;48,-54;46,-55;44,-57;42,-58;39,-60;37,-61;35,-63;33,-64;30,-66;27,-67;24,-68;21,-69;18,-70;16,-70;14,-71;10,-71;7,-72;5,-72;2,-72;0,-72;-2,-72;-5,-72;-7,-72;-10,-72;-12,-72;-15,-71;-18,-71;-20,-70;-23,-69;-25,-69;-28,-68;-31,-67;-33,-66;-35,-65;-37,-63;-39,-62;-42,-61;-44,-59;-47,-58;-48,-57;-51,-55;-52,-54;
% ];

N = size(shift_ex, 1);

% x_pos = round( x_pos * N );
% if x_pos == N
%     x_pos = 0;
% end
% x_pos = x_pos + 1;
% the_shift_ex = shift_ex(x_pos, 1:2);

x_pos_N   = x_pos * 180;
x_pos_N_l = floor(x_pos_N);
x_pos_N_r = ceil (x_pos_N);
if x_pos_N_l==x_pos_N_r
    % 直接给出结果 give the result directly
    x_idx = x_pos_N_l + 1;
    the_shift_ex = shift_ex(x_idx, 1:2);
else
    % 简单的线性插值 simple linear interpolation
    x_idx_l = x_pos_N_l + 1;
    x_idx_r = x_pos_N_r + 1; 
    if x_idx_r>N
        x_idx_r = x_idx_r - N;
    end
    
    the_shift_ex_l = shift_ex(x_idx_l, 1:2);
    the_shift_ex_r = shift_ex(x_idx_r, 1:2);
    dist_l = x_pos_N   - x_pos_N_l;
    dist_r = x_pos_N_r - x_pos_N  ; 
    
    the_shift_ex = the_shift_ex_r * dist_l + the_shift_ex_l * dist_r;
end

the_shift_ex = round(the_shift_ex);

end

