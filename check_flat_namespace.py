#!/usr/bin/env python3
"""Check if a Mach-O binary has flat namespace set"""
import struct
import sys

def check_macho_flags(filename):
    """Check Mach-O flags"""
    try:
        with open(filename, 'rb') as f:
            magic = struct.unpack('<I', f.read(4))[0]
            
            if magic == 0xfeedfacf:  # MH_MAGIC_64
                f.read(16)  # Skip cputype, cpusubtype, filetype, ncmds, sizeofcmds
                flags = struct.unpack('<I', f.read(4))[0]
                
                MH_TWOLEVEL = 0x80
                MH_FORCE_FLAT = 0x100
                
                has_flat = bool(flags & MH_FORCE_FLAT)
                has_twolevel = bool(flags & MH_TWOLEVEL)
                
                if has_flat:
                    print(f"✅ {filename.split('/')[-1]}: FLAT NAMESPACE (0x{flags:08x})")
                    return True
                elif has_twolevel:
                    print(f"❌ {filename.split('/')[-1]}: TWO-LEVEL namespace (0x{flags:08x})")
                    return False
                else:
                    print(f"⚠️  {filename.split('/')[-1]}: Neither set explicitly (0x{flags:08x})")
                    return False
            else:
                print(f"❌ Not a Mach-O 64-bit file")
                return False
    except Exception as e:
        print(f"❌ Error reading {filename}: {e}")
        return False

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: check_flat_namespace.py <file.so> [file2.so ...]")
        sys.exit(1)
    
    results = []
    for filename in sys.argv[1:]:
        results.append(check_macho_flags(filename))
    
    if all(results):
        print("\n✅ All files have flat namespace!")
        sys.exit(0)
    else:
        print("\n❌ Some files missing flat namespace!")
        sys.exit(1)
