---
title: Uutolsó eredményeim
---

# utolsó eredményeim

DDD szimulációk optimalizálása és korrelációk vizsgálata

- [utolsó eredményeim](#utolsó-eredményeim)
- [Fordítási optimalizáció](#fordítási-optimalizáció)
  - [Mit csinálnak ezek a kapcsolók?](#mit-csinálnak-ezek-a-kapcsolók)
- [Feszültségtér optimalizációja](#feszültségtér-optimalizációja)
- [Korrelációk](#korrelációk)

# Fordítási optimalizáció

Eredeti Makefileon további beállítások:

```
set(USR_COMP_OPTIONS "-Wall -pedantic -Wextra -Wno-unknown-pragmas -Ofast -flto  -ffat-lto-objects -march=native -ffast-math")
```

## Mit csinálnak ezek a kapcsolók?

- ha bármi gyanús van, akkor álljon le:

  `-Wall -pedantic -Wextra`
- hadd használjam a commenteket region-ok definiálására (code structure), és lehessenek OS specifikus beállítások:
  
  `-Wno-unknown-pragmas`

- `-Ofast`: optimalizáld futási időre (cserébe lehet, hogy nagyobb lesz, vagy sokáig fordul)
- `flto -ffat-lto-objects`: link time optimalizáció, szélsőséges esetben akár 2-es szorzót is hozhat. A második elvileg nem kell, de mintha mégis hasznos lett volna.
- `-march=native`: ha CPU-specifikus gyorsítást lehet elérni, tedd meg (cserébe nem fog más CPU-n futni)
- `-ffast-math`: szegd meg az ieee754 szabványt a gyorsabb futásért cserébe
  - 0-val való osztás, inf és nan-ok más kezelése
  - átzárójelezés: `(a * b) * c` helyett `a * (b * c)`. Hasznos pl: `(1 * a) * 1/a` helyett `1 * (a * 1/a) = 1`
  - olcsóbb függvény használata: `pow(2.5,3)` helyett `2.5 * 2.5 * 2.5`-öt használ.

# Feszültségtér optimalizációja



# Korrelációk

