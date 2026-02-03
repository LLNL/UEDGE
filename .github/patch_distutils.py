"""
Patch distutils to remove -undefined dynamic_lookup on macOS and ensure flat_namespace.
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
        """Remove -undefined dynamic_lookup and add flat_namespace"""
        _original_customize_compiler(compiler)
        
        if hasattr(compiler, 'linker_so'):
            # Remove -undefined dynamic_lookup
            if '-undefined' in compiler.linker_so:
                idx = compiler.linker_so.index('-undefined')
                # Remove both -undefined and dynamic_lookup
                if idx + 1 < len(compiler.linker_so):
                    compiler.linker_so.pop(idx)  # Remove -undefined
                    compiler.linker_so.pop(idx)  # Remove dynamic_lookup
                    print(f"  Removed -undefined dynamic_lookup from linker_so")
            
            # Ensure flat_namespace is present (with -Wl, prefix)
            if '-Wl,-flat_namespace' not in compiler.linker_so:
                compiler.linker_so.extend(['-Wl,-flat_namespace', '-Wl,-U,___emutls_get_address'])
                print(f"  Added -Wl,-flat_namespace to linker_so")
            
            print(f"  Final linker_so: {' '.join(compiler.linker_so)}")
        
        return compiler
    
    sysconfig.customize_compiler = custom_customize_compiler
    distutils.sysconfig.customize_compiler = custom_customize_compiler
    
    # Also patch the LDSHARED config var
    original_get_config_vars = sysconfig.get_config_vars
    
    def patched_get_config_vars(*args):
        """Patch LDSHARED to remove -undefined dynamic_lookup"""
        result = original_get_config_vars(*args)
        if isinstance(result, dict) and 'LDSHARED' in result:
            ldshared = result['LDSHARED']
            if '-undefined dynamic_lookup' in ldshared:
                ldshared = ldshared.replace('-undefined dynamic_lookup', '')
                ldshared += ' -Wl,-flat_namespace -Wl,-U,___emutls_get_address'
                result['LDSHARED'] = ldshared
                print(f"  Patched LDSHARED: {ldshared}")
        return result
    
    sysconfig.get_config_vars = patched_get_config_vars
    distutils.sysconfig.get_config_vars = patched_get_config_vars
    
    print("=== Distutils patching complete ===")
