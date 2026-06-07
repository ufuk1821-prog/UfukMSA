# UfukMSA

MAFFT tabanlı Çoklu Dizi Hizalama (Multiple Sequence Alignment) kütüphanesi.

## Kurulum

```
pip install UfukMSA
```

## Kullanım

```python
from ufukmsa import mafft

sequences = {
    "seq1": "ACGTACGT",
    "seq2": "ACGTTCGT",
    "seq3": "ACGTACG",
}

result = mafft(sequences)
for name, aligned in result.items():
    print(f"{name}: {aligned}")
```

## Algoritma

Bu kütüphane MAFFT algoritmasının basitleştirilmiş bir uygulamasıdır:
1. Tüm dizi çiftleri arasında Needleman-Wunsch hizalaması yapılır
2. Benzerlik skorlarından uzaklık matrisi oluşturulur
3. UPGMA yöntemiyle guide tree (rehber ağacı) oluşturulur
4. Guide tree sırasına göre profil-profil hizalaması yapılır
