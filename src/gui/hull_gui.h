/**
 * convex hull test
 * @author Tobias Weber
 * @date apr-2021
 * @license see 'LICENSE' file
 */

#ifndef __GEO_HULL_H__
#define __GEO_HULL_H__

#include <QDialog>
#include <QTableWidget>
#include <QToolButton>
#include <QPlainTextEdit>
#include <QMenu>

#include <vector>
#include <sstream>
#include <iostream>



template<class T = double>
class NumericTableWidgetItem : public QTableWidgetItem
{
public:
	NumericTableWidgetItem(const T& val)
		: QTableWidgetItem(std::to_string(val).c_str()),
		m_val{val}
	{}

	virtual ~NumericTableWidgetItem() = default;


	virtual bool operator<(const QTableWidgetItem& item) const override
	{
		std::istringstream istr1{this->text().toStdString()};
		std::istringstream istr2{item.text().toStdString()};

		T val1{}, val2{};
		istr1 >> val1;
		istr2 >> val2;

		return val1 < val2;
	}


	virtual void setData(int itemdatarole, const QVariant& var) override
	{
		if(itemdatarole == Qt::EditRole)
			m_val = var.value<T>();
		QTableWidgetItem::setData(itemdatarole, std::to_string(m_val).c_str());
	}


	virtual QTableWidgetItem* clone() const override
	{
		return new NumericTableWidgetItem<T>(m_val);
	}


	T GetValue() const { return m_val; }


private:
	T m_val{};
};



class HullDlg : public QDialog
{
public:
	HullDlg(QWidget* pParent = nullptr);
	HullDlg(const HullDlg&) = default;

	virtual ~HullDlg() = default;

	void CalculateHull();


protected:
	QWidget *m_pTabWidget{};
	QTableWidget *m_pTab{};
	QPlainTextEdit *m_pEdit{};

	QToolButton *m_pTabBtnAdd{};
	QToolButton *m_pTabBtnDel{};
	QToolButton *m_pTabBtnUp{};
	QToolButton *m_pTabBtnDown{};

	QMenu *m_pTabContextMenu{};


protected:
	void AddTabItem(int row = -1);
	void DelTabItem();
	void MoveTabItemUp();
	void MoveTabItemDown();

	void TableCellChanged(int rowNew, int colNew, int rowOld, int colOld);
	void ShowTableContextMenu(const QPoint& pt);


private:
	std::vector<int> GetSelectedRows(bool sort_reversed = false) const;


private:
	int m_iCursorRow = -1;

};


#endif
