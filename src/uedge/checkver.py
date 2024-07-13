pkg = "uedge"

try:
    import json
    import urllib.request

    try:
        import importlib.metadata

        thisver = importlib.metadata.version(pkg)
    except:
        import pkg_resources

        thisver = pkg_resources.get_distribution(pkg).version

    contents = urllib.request.urlopen("https://pypi.org/pypi/" + pkg + "/json").read()
    data = json.loads(contents.decode())
    thatver = data["info"]["version"]

    if thisver < thatver:
        print()
        print("Uedge version " + thisver + ", an update is available to " + thatver)
        print()

except Exception as err:
    print()
    print("Error checking pypi version: {}".format(err))
    print()
