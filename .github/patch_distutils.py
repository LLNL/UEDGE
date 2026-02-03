"""
Patch distutils to ensure flat_namespace on macOS for OpenMP emutls compatibility.
This should be imported at the very top of setup.py before anything else.
"""
import sys
import os

if sys.platform == 'darwin':
    print("=== Patching distutils for macOS flat namespace ===")
    
    # Patch distutils/sysconfig
    import distutils.sysconfig
    from distutils import sysconfig
    
    _original_customize_compiler = sysconfig.customize_compiler
    
    def custom_customize_compiler(compiler):
        """Ensure flat_namespace is set for linking (keep -undefined dynamic_lookup)"""
        _original_customize_compiler(compiler)
        
        if hasattr(compiler, 'linker_so'):
            # Add flat_namespace AFTER -undefined dynamic_lookup (order matters!)
            # This ensures flat namespace takes effect
            flags_to_add = ['-Wl,-flat_namespace', '-Wl,-U,___emutls_get_address']
            
            for flag in flags_to_add:
                if flag not in compiler.linker_so:
                    compiler.linker_so.append(flag)
            
            print(f"  Patched linker_so: {' '.join(compiler.linker_so)}")
        
        # IMPORTANT: Don't add linker flags to compiler_so (compilation, not linking)
        # Only linker_so should have these flags
        
        return compiler
    
    sysconfig.customize_compiler = custom_customize_compiler
    distutils.sysconfig.customize_compiler = custom_customize_compiler
    
    print("=== Distutils patching complete ===")
