classdef Mmaptx < handle
    % Matlab structure for memmap transfering with python
    %   Class that enables memory mapped files data transfer method with
    %   the respective python library.
    %
    % WARNING: This was slower than direct implementation. It's better to
    % use the memory map directly instead through this interface if speed
    % is important.

    properties
        mm_in
        mm_out
        blocking
        sync = 0;
    end

    methods
        function obj = Mmaptx(name, ftype, blocking)
            %Mmaptx Construct an instance of this class
            %   Initializes the memory mapped files for transfer and the
            %   class parameters/properties

            obj.mm_in = memmapfile(strcat(name,'_mmap_out.dat'),'Format', ftype);

            obj.mm_out = memmapfile(strcat(name,'_mmap_in.dat'),'Format', ftype);
            obj.mm_out.Writable = true;
            obj.blocking = blocking;
        end

        function write(obj, input)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.mm_out.Data(2:end) = input;
            obj.mm_out.Data(1) = obj.sync;
            
        end

        function data = read(obj)
            %read Reads data from the memory mapped file
            %   Read data from the memory mapped file
            %   If blocking is enabled, it will read only when the sync
            %   flag changes in the input file.
            %   Returns the read data and updates the sync flag for sending

            if obj.blocking == true
                while obj.mm_in.Data(1) == obj.sync; end % blocking wait
            end
            obj.sync = obj.mm_in.Data(1);
            data = obj.mm_in.Data(2:end);
        end
    end
end