import mmap
import struct
import atexit
import os

# Library that can communicate with other programs through memory mapped files
# It writes to an output memory mapped file and it reads from a different input memory mapped file
# It can also sync data with the other program via a sync variable that is passed to the other
# program and input data will only be read unless that flag is sent back in the input memory map file (if blocking read is enabled)
class Mmaptx:

    def __init__(self, name = None, format_type = "d", in_size = 1, out_size = 1, blocking = True):

        atexit.register(self.cleanup)

        # Save configuration data
        self.ftype = format_type
        self.ftype_size = struct.calcsize(format_type)
        self.blocking = blocking
        self.in_size = in_size
        self.out_size = out_size

        if name == None:
            exit("Name must be specified.")
        
        self.filename = name
        self.filename_in = name + "_mmap_in.dat"
        self.filename_out = name + "_mmap_out.dat"

        # Create or clean files that will be used as memory map input and output
        f = open(self.filename_in, "w+b")
        f.write(struct.pack(self.ftype * (self.in_size+1), *([0]*(self.in_size+1))))
        f.flush()
        f.close()

        f = open(self.filename_out, "w+b")
        f.write(struct.pack(self.ftype * (self.out_size+1), *([0]*(self.out_size+1))))
        f.flush()
        f.close()

        del f

        # Crete the memory maps
        self.f_in = open(self.filename_in, "r+b")
        self.mm_in = mmap.mmap(self.f_in.fileno(), 0)

        self.f_out = open(self.filename_out, "a+b")
        self.mm_out = mmap.mmap(self.f_out.fileno(), 0)

        self.sync = 0.0
    
    # Write to output memmap
    def write(self, *data):
        
        out_size_l = len(data)
        self.mm_out.seek(self.ftype_size)
        self.mm_out.write(struct.pack(self.ftype * out_size_l, *data))
        
        self.sync = self.sync
        self.mm_out.seek(0)
        self.mm_out.write(struct.pack(self.ftype, self.sync))

    # Blocking read from input memmap
        # Waits until the other program has set the sync flag to the same value sent
        # This means data is ready for reading
    def read(self):

        self.mm_in.seek(0)
        sync_in = struct.unpack(self.ftype, self.mm_in.read(self.ftype_size))

        if self.blocking:
            while self.sync == sync_in[0]:
                self.mm_in.seek(0)
                sync_in = struct.unpack(self.ftype, self.mm_in.read(self.ftype_size))
            self.sync = sync_in[0]

        data_in = struct.unpack(self.ftype * self.in_size, self.mm_in.read(self.ftype_size * self.in_size))
            
        return data_in
    
    # Close files and memmaps
    def cleanup(self):

        self.mm_out.close()
        self.f_out.close()
        self.mm_in.close()
        self.f_in.close()
