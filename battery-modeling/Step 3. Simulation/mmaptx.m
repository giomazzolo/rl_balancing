classdef mmaptx
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        mm_in;
        mm_out;
    end

    methods
        function obj = mmaptx(name, ftype)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            this.m_in = memmapfile(strcat(name,'simCell_mmap_out.dat'),'Format', ftype);

            this.m_out = memmapfile(strcat(name,'_mmap_in.dat'),'Format', ftype);
            this.m_out.Writable = true;
        end

        function write(inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            this.m_out.Data(3) = inputArg.vk;
            this.m_out.Data(2) = inputArg.state;
            this.m_out.Data(1) = inputArg.sync;
        end
    end
end