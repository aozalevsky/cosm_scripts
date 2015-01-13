part = ['H', 'PT', 'T', 'T1', 'T2', 'T3', 'T4', 'T5',
        'T6', 'B', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
        'N', 'S']
dih = []

for i, a in enumerate(part):
    b = a + 'A'
    part[i] = b

for i1 in part:
    for i2 in part:
        for i3 in part:
            for i4 in part:
                a = [i1, i2, i3, i4]
                if a not in dih and a[::-1] not in dih:
                    dih.append(a)
for a in dih:
    print a[0], a[1], a[2], a[3], 1, 0.0, 1.00, 1
