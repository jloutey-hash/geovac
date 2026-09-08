import subprocess, sys
for n in (7, 8, 9, 10, 11, 12, 13, 14):
    subprocess.run([sys.executable, "-c",
                    f"import sys;sys.path.insert(0,'.');sys.argv=['x','{n}','24000'];"
                    "exec(open('debug/p60_boxrule.py').read())"], check=True)
