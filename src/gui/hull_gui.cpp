/**
 * convex hull test
 * @author Tobias Weber
 * @date apr-2021
 * @license: see 'LICENSE' file
 */

#include "hull_gui.h"

#include <QApplication>
#include <QGridLayout>
#include <QHeaderView>
#include <QSpinBox>

#include <locale>
#include <iostream>

#include "../src/geo_algos.h"

using t_real = double;
using t_vec = m::vec<t_real, std::vector>;
using t_mat = m::mat<t_real, std::vector>;



// ----------------------------------------------------------------------------
HullDlg::HullDlg(QWidget* pParent) : QDialog{pParent}
{
	setWindowTitle("Convex Hull");
	m_pTabWidget = new QWidget(this);


	// table
	m_pTab = new QTableWidget(m_pTabWidget);
	m_pTab->setShowGrid(true);
	m_pTab->setSortingEnabled(true);
	m_pTab->setMouseTracking(true);
	m_pTab->setSelectionBehavior(QTableWidget::SelectRows);
	m_pTab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_pTab->setContextMenuPolicy(Qt::CustomContextMenu);

	m_pTab->setColumnCount(3);
	m_pTab->setHorizontalHeaderItem(0, new QTableWidgetItem{"x"});
	m_pTab->setHorizontalHeaderItem(1, new QTableWidgetItem{"y"});
	m_pTab->setHorizontalHeaderItem(2, new QTableWidgetItem{"z"});

	m_pTab->horizontalHeader()->setDefaultSectionSize(200);
	m_pTab->verticalHeader()->setDefaultSectionSize(32);
	m_pTab->verticalHeader()->setVisible(false);

	m_pTab->setColumnWidth(0, 150);
	m_pTab->setColumnWidth(1, 150);
	m_pTab->setColumnWidth(2, 150);


	// text edit
	m_pEdit = new QPlainTextEdit(m_pTabWidget);
	m_pEdit->setReadOnly(true);


	// buttons
	m_pTabBtnAdd = new QToolButton(m_pTabWidget);
	m_pTabBtnDel = new QToolButton(m_pTabWidget);
	m_pTabBtnUp = new QToolButton(m_pTabWidget);
	m_pTabBtnDown = new QToolButton(m_pTabWidget);

	m_pTabBtnAdd->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});
	m_pTabBtnDel->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});
	m_pTabBtnUp->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});
	m_pTabBtnDown->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});

	m_pTabBtnAdd->setText("\xe2\x8a\x95");
	m_pTabBtnDel->setText("\xe2\x8a\x96");
	m_pTabBtnUp->setText("\342\206\221");
	m_pTabBtnDown->setText("\342\206\223");

	m_pTabBtnAdd->setToolTip("Add item.");
	m_pTabBtnDel->setToolTip("Delete selected item(s).");
	m_pTabBtnUp->setToolTip("Move selected item(s) up.");
	m_pTabBtnDown->setToolTip("Move selected item(s) down.");


	// table grid
	int y = 0;
	auto pTabGrid = new QGridLayout(m_pTabWidget);
	pTabGrid->setSpacing(2);
	pTabGrid->setContentsMargins(4,4,4,4);
	pTabGrid->addWidget(m_pTab, y++,0,1,5);
	pTabGrid->addWidget(m_pEdit, y++,0,1,5);
	pTabGrid->addWidget(m_pTabBtnAdd, y,0,1,1);
	pTabGrid->addWidget(m_pTabBtnDel, y,1,1,1);
	pTabGrid->addItem(new QSpacerItem(4, 4, QSizePolicy::Expanding, QSizePolicy::Minimum), y,2,1,1);
	pTabGrid->addWidget(m_pTabBtnUp, y,3,1,1);
	pTabGrid->addWidget(m_pTabBtnDown, y++,4,1,1);


	// table context CustomContextMenu
	m_pTabContextMenu = new QMenu(m_pTab);
	m_pTabContextMenu->addAction("Add Item Before", this, [this]() { this->AddTabItem(-2); });
	m_pTabContextMenu->addAction("Add Item After", this, [this]() { this->AddTabItem(-3); });
	m_pTabContextMenu->addAction("Delete Item", this, &HullDlg::DelTabItem);


	// signals
	connect(m_pTabBtnAdd, &QToolButton::clicked, this, [this]() { this->AddTabItem(-1); });
	connect(m_pTabBtnDel, &QToolButton::clicked, this, &HullDlg::DelTabItem);
	connect(m_pTabBtnUp, &QToolButton::clicked, this, &HullDlg::MoveTabItemUp);
	connect(m_pTabBtnDown, &QToolButton::clicked, this, &HullDlg::MoveTabItemDown);
	connect(m_pTab, &QTableWidget::currentCellChanged, this, &HullDlg::TableCellChanged);
	connect(m_pTab, &QTableWidget::customContextMenuRequested, this, &HullDlg::ShowTableContextMenu);


	// main grid
	auto pDlgGrid = new QGridLayout(this);
	pDlgGrid->setSpacing(2);
	pDlgGrid->setContentsMargins(4,4,4,4);
	pDlgGrid->addWidget(m_pTabWidget, 0,0,1,1);
}


void HullDlg::CalculateHull()
{
	using namespace m_ops;

	try
	{
		m_pEdit->clear();

		int dim = m_pTab->columnCount();
		int rows = m_pTab->rowCount();

		if(rows < 4)
			return;

		// get vertices
		std::vector<t_vec> vertices;

		for(int row=0; row<rows; ++row)
		{
			t_vec vertex = m::create<t_vec>(dim);

			for(int col=0; col<dim; ++col)
			{
				auto widget = dynamic_cast<NumericTableWidgetItem<t_real>*>(m_pTab->item(row, col));
				vertex[col] = widget->GetValue();
			}

			vertices.emplace_back(std::move(vertex));
		}

		// calculate hull
		std::vector<std::vector<t_vec>> hull;
		std::tie(std::ignore, hull, std::ignore) = g::calc_delaunay<t_vec>(dim, vertices, true);


		// output results
		std::ostringstream ostr;

		for(std::size_t vertidx=0; vertidx<vertices.size(); ++vertidx)
		{
			ostr << "Vertex " << (vertidx+1) << ": "
				<< vertices[vertidx] << "\n";
		}

		ostr << "\n";

		for(std::size_t polyidx=0; polyidx<hull.size(); ++polyidx)
		{
			const auto& poly = hull[polyidx];

			ostr << "Polygon " << (polyidx+1) << ":\n";

			for(const auto& vertex : poly)
				ostr << "\t" << vertex << "\n";

			ostr << "\n";
		}

		m_pEdit->setPlainText(ostr.str().c_str());
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}
}


void HullDlg::AddTabItem(int row)
{
	// append to end of table
	if(row == -1)
		row = m_pTab->rowCount();

	// use row from member variable
	else if(row == -2 && m_iCursorRow >= 0)
		row = m_iCursorRow;

	// use row from member variable +1
	else if(row == -3 && m_iCursorRow >= 0)
		row = m_iCursorRow + 1;

	m_pTab->setSortingEnabled(false);
	m_pTab->insertRow(row);

	m_pTab->setItem(row, 0, new NumericTableWidgetItem<t_real>(0));
	m_pTab->setItem(row, 1, new NumericTableWidgetItem<t_real>(0));
	m_pTab->setItem(row, 2, new NumericTableWidgetItem<t_real>(0));

	m_pTab->scrollToItem(m_pTab->item(row, 0));
	m_pTab->setCurrentCell(row, 0);
	m_pTab->setSortingEnabled(true);
}


void HullDlg::DelTabItem()
{
	// if nothing is selected, clear all items
	if(m_pTab->selectedItems().count() == 0)
	{
		m_pTab->clearContents();
		m_pTab->setRowCount(0);
	}


	for(int row : GetSelectedRows(true))
		m_pTab->removeRow(row);
}


void HullDlg::MoveTabItemUp()
{
	m_pTab->setSortingEnabled(false);

	auto selected = GetSelectedRows(false);
	for(int row : selected)
	{
		if(row == 0)
			continue;

		auto *item = m_pTab->item(row, 0);
		if(!item || !item->isSelected())
			continue;

		m_pTab->insertRow(row-1);
		for(int col=0; col<m_pTab->columnCount(); ++col)
			m_pTab->setItem(row-1, col, m_pTab->item(row+1, col)->clone());
		m_pTab->removeRow(row+1);
	}

	for(int row=0; row<m_pTab->rowCount(); ++row)
	{
		if(auto *item = m_pTab->item(row, 0);
			item && std::find(selected.begin(), selected.end(), row+1) != selected.end())
		{
			for(int col=0; col<m_pTab->columnCount(); ++col)
				m_pTab->item(row, col)->setSelected(true);
		}
	}
}


void HullDlg::MoveTabItemDown()
{
	m_pTab->setSortingEnabled(false);

	auto selected = GetSelectedRows(true);
	for(int row : selected)
	{
		if(row == m_pTab->rowCount()-1)
			continue;

		auto *item = m_pTab->item(row, 0);
		if(!item || !item->isSelected())
			continue;

		m_pTab->insertRow(row+2);
		for(int col=0; col<m_pTab->columnCount(); ++col)
			m_pTab->setItem(row+2, col, m_pTab->item(row, col)->clone());
		m_pTab->removeRow(row);
	}

	for(int row=0; row<m_pTab->rowCount(); ++row)
	{
		if(auto *item = m_pTab->item(row, 0);
			item && std::find(selected.begin(), selected.end(), row-1) != selected.end())
		{
			for(int col=0; col<m_pTab->columnCount(); ++col)
				m_pTab->item(row, col)->setSelected(true);
		}
	}
}


std::vector<int> HullDlg::GetSelectedRows(bool sort_reversed) const
{
	std::vector<int> vec;
	vec.reserve(m_pTab->selectedItems().size());

	for(int row=0; row<m_pTab->rowCount(); ++row)
	{
		if(auto *item = m_pTab->item(row, 0); item && item->isSelected())
			vec.push_back(row);
	}

	if(sort_reversed)
	{
		std::stable_sort(vec.begin(), vec.end(), [](int row1, int row2)
		{ return row1 > row2; });
	}

	return vec;
}


void HullDlg::TableCellChanged(int rowNew, int colNew, int rowOld, int colOld)
{
	CalculateHull();
}


void HullDlg::ShowTableContextMenu(const QPoint& pt)
{
	const auto* item = m_pTab->itemAt(pt);
	if(!item)
		return;

	m_iCursorRow = item->row();
	auto ptGlob = m_pTab->mapToGlobal(pt);
	ptGlob.setY(ptGlob.y() + m_pTabContextMenu->sizeHint().height()/2);
	m_pTabContextMenu->popup(ptGlob);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
int main(int argc, char** argv)
{
	auto app = std::make_unique<QApplication>(argc, argv);
	std::ios_base::sync_with_stdio(false);
	::setlocale(LC_ALL, "C");
	std::locale::global(std::locale("C"));
	QLocale::setDefault(QLocale::C);

	auto dlg = std::make_unique<HullDlg>(nullptr);
	dlg->resize(600, 500);
	dlg->show();

	return app->exec();
}
// ----------------------------------------------------------------------------
