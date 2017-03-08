__author__ = 'anna'

if __name__ == '__main__':
    import datetime
    import platform
    machine = platform.system() + '_' + platform.machine()
    now = str(datetime.datetime.now())[:10]
    with open('head.py', 'r') as rp:
        content = rp.readlines()
    for j, line in enumerate(content):
        if line.startswith('    #                  by Anna V. Luebben ('):
            old = line.partition('(')[2].partition(')')[0]
            now = machine + '@' + now
            if old == now:
                    now += '#2'
                    break
            # if '_' in old:
            #     i = int(old.partition('#')[2]) + 1
            #     now = now + '_#{}'.format(i)

            break
    now += ')'
    line = '    #                  by Anna V. Luebben ({: <32}#\n'.format(now)
    content[j] = line
    with open('head.py', 'w') as wp:
        for line in content:
            wp.write(line)

