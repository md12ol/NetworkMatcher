import os

inp = "../../Conferences and Papers/2023 CIBCB/WSDA/BitOutCheck/"


def combine(dir: str, start: str, end: str):
    if os.path.exists(dir + start + end):
        os.remove(dir + start + end)
        pass
    out_file = open(dir + start + end, "a")
    for idx in range(1,31):
        with open(dir + start + str(idx).zfill(2) + end, "r") as f:
            out_file.write(f.read())
            pass
        pass
    out_file.close()
    pass 


def main():
    folder_names = os.listdir(inp)
    for fold in folder_names:
        combine(inp + fold + "/", "best", ".lint")
        # combine(inp + fold + "/", "patient0", ".dat")
        pass
    pass


main()