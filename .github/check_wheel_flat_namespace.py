#!/usr/bin/env python3
"""
Check if wheel modules have flat namespace set.
Usage: python check_wheel_flat_namespace.py <wheel_file>
"""
import struct
import glob
import sys
import zipfile
import tempfile
import os
import shutil

def check_flat_namespace(so_file):
    """Check if a .so file has flat namespace flag set."""
    try:
        with open(so_file, 'rb') as f:
            magic = struct.unpack('<I', f.read(4))[0]
            if magic == 0xfeedfacf:  # MH_MAGIC_64
                f.read(16)  # Skip to flags
                flags = struct.unpack('<I', f.read(4))[0]
                
                MH_FORCE_FLAT = 0x100
                has_flat = bool(flags & MH_FORCE_FLAT)
                
                return has_flat, flags
    except Exception as e:
        print(f"  Error reading {so_file}: {e}")
        return None, None
    
    return None, None

def main():
    if len(sys.argv) < 2:
        print("Usage: python check_wheel_flat_namespace.py <wheel_file>")
        sys.exit(1)
    
    wheel_file = sys.argv[1]
    
    if not os.path.exists(wheel_file):
        print(f"Error: Wheel file not found: {wheel_file}")
        sys.exit(1)
    
    print(f"Checking wheel: {wheel_file}")
    print("=" * 60)
    
    # Extract wheel to temp directory
    tmpdir = tempfile.mkdtemp()
    
    try:
        with zipfile.ZipFile(wheel_file, 'r') as z:
            z.extractall(tmpdir)
        
        # Check critical modules
        modules_to_check = ['ppppy', 'bbbpy', 'apipy']
        results = []
        
        for module_name in modules_to_check:
            pattern = os.path.join(tmpdir, 'uedge', f'*{module_name}*.so')
            so_files = glob.glob(pattern)
            
            if not so_files:
                print(f"⚠️  {module_name}: Not found")
                continue
            
            for so_file in so_files:
                has_flat, flags = check_flat_namespace(so_file)
                
                if has_flat is None:
                    print(f"❌ {os.path.basename(so_file)}: Not a Mach-O file or error")
                    results.append(False)
                elif has_flat:
                    print(f"✅ {os.path.basename(so_file)}: FLAT NAMESPACE (0x{flags:08x})")
                    results.append(True)
                else:
                    print(f"❌ {os.path.basename(so_file)}: TWO-LEVEL namespace (0x{flags:08x})")
                    results.append(False)
        
        print("=" * 60)
        
        if not results:
            print("⚠️  No modules found to check")
            sys.exit(0)
        elif all(results):
            print("✅ SUCCESS: All checked modules have flat namespace!")
            sys.exit(0)
        else:
            print("❌ FAILURE: Some modules missing flat namespace!")
            sys.exit(1)
            
    finally:
        # Cleanup
        shutil.rmtree(tmpdir, ignore_errors=True)

if __name__ == '__main__':
    main()
