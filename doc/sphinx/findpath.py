import os

def FindPath(relpath):
    def Try(abspath, *path):
        s = os.path.abspath(os.path.join(*path))
        if os.path.exists(s):
            assert not abspath, \
                "Ambiguous path: target found in '{:}' and '{:}'".format(
                        abspath, s)
            return s
        return abspath

    r = relpath
    repo = os.path.abspath(os.path.join(os.getcwd(), "../../"))
    a = "" # absolute path
    a = Try(a, repo, r)
    a = Try(a, repo, "deploy", r)
    a = Try(a, repo, "deploy/lib", r)
    a = Try(a, repo, "deploy/tool", r)
    a = Try(a, repo, "deploy/wrap", r)
    a = Try(a, repo, "src", r)
    return a
