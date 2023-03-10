def update_progress(tot):
    # Could be shared with extract features?
    prog = 0
    while True:
        yield prog
        prog += 1
        if prog % 20 == 0: print("\t Progress: {}/{}".format(prog, tot))
