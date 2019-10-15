from PyQt5.QtWidgets import QWidget, QApplication, QComboBox, QFormLayout

class Widget(QWidget):
    def __init__(self, parent = None):
        super().__init__(parent)
        self.categories = {'animals':['cat', 'dog', 'parrot', 'fish'],
                           'flowers':['daisies', 'tulips', 'daffodils', 'roses'],
                           'colors':['red', 'orange', 'blue', 'purple']}
        self.cat_combobox = QComboBox(self)
        self.item_combobox = QComboBox(self)

        self.cat_combobox.setEditable(False)
        self.item_combobox.setEditable(False)

        self.cat_combobox.currentTextChanged.connect(self.set_category)
        self.cat_combobox.addItems(sorted(self.categories.keys()))

        form_layout = QFormLayout(self)
        form_layout.addRow('Category', self.cat_combobox)
        form_layout.addRow('Items', self.item_combobox)

    def set_category(self, text):
        self.item_combobox.clear()
        self.item_combobox.addItems(self.categories.get(text, []))

if __name__ == '__main__':
    import sys

    app = QApplication(sys.argv)
    window = Widget()
    window.show()
    # app.exec()