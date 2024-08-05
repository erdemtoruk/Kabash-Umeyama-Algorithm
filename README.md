# Kabash-Umeyama-Algorithm

Bu proje, Python kullanarak Kabash-Umeyama algoritmasını gerçekleştirmektedir. Python.linalg içerisindeki fonksiyonlar kullanılmadan tekrar yazılmıştır.

## Algoritmanın Amacı

Kabash-Umeyama algoritması, bir kaynak nokta kümesini bir hedef nokta kümesiyle hizalamak için kullanılır. Bu işlem, iki nokta kümesi arasındaki en uygun rotasyon, çeviri ve ölçekleme parametrelerini hesaplamayı içerir. Algoritma, her iki nokta kümesinin aynı boyutta ve sıralı olduğunu varsayar.

### Algoritmanın Adımları

1. **Dosyaları okuma**: Nokta kümelerini içeren dosyalar okunur.
2. **Ortalama Hesaplama**: Her iki nokta kümesi için ortalama noktalar hesaplanır.
3. **Kovaryans Matrisi**: Kaynak ve hedef noktalar arasındaki kovaryans matrisi hesaplanır.
4. **Singular Value Decomposition (SVD)**: Kovaryans matrisi üzerinde SVD uygulanır.
5. **Rotasyon Matrisi**: SVD sonuçlarından optimal rotasyon matrisi hesaplanır.
6. **Çeviri Vektörü**: Çeviri vektörü hesaplanır.
7. **Hizalama**: Kaynak noktalar, hesaplanan rotasyon, ölçek ve çeviri parametreleri kullanılarak hedef noktalara hizalanır.

#### Lisans

Bu proje MIT lisansı ile lisanslanmıştır.