/**
 * vis test program
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 11-Nov-2020
 * @license: see 'LICENSE' file
 */

#include "vis_gui.h"

#include <QApplication>
#include <QMenuBar>
#include <QLabel>
#include <QStatusBar>
#include <QMouseEvent>
#include <QFileDialog>
#include <QSvgGenerator>
#include <QMessageBox>
#include <QSettings>

#include <locale>
#include <memory>
#include <array>
#include <vector>
#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace ptree = boost::property_tree;



// ----------------------------------------------------------------------------

Vertex::Vertex(const QPointF& pos, double rad) : m_rad{rad}
{
	setPos(pos);
	setFlags(flags() | QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable);
}


Vertex::~Vertex()
{
}


QRectF Vertex::boundingRect() const
{
	return QRectF{-m_rad/2., -m_rad/2., m_rad, m_rad};
}


void Vertex::paint(QPainter* painter, const QStyleOptionGraphicsItem*, QWidget*)
{
	std::array<QColor, 2> colours =
	{
		QColor::fromRgbF(0.,0.,1.),
		QColor::fromRgbF(0.,0.,0.),
	};

	QRadialGradient grad{};
	grad.setCenter(0., 0.);
	grad.setRadius(m_rad);

	for(std::size_t col=0; col<colours.size(); ++col)
		grad.setColorAt(col/double(colours.size()-1), colours[col]);

	painter->setBrush(grad);
	painter->setPen(*colours.rbegin());

	painter->drawEllipse(-m_rad/2., -m_rad/2., m_rad, m_rad);
}

// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------

VisView::VisView(QGraphicsScene *scene, QWidget *parent) : QGraphicsView(scene, parent),
	m_scene{scene}
{
	setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
	setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

	setInteractive(true);
	setMouseTracking(true);

	//scale(1., -1.);

	setBackgroundBrush(QBrush{QColor::fromRgbF(0.95, 0.95, 0.95, 1.)});
}


VisView::~VisView()
{
}


void VisView::resizeEvent(QResizeEvent *evt)
{
	QPointF pt1{mapToScene(QPoint{0,0})};
	QPointF pt2{mapToScene(QPoint{evt->size().width(), evt->size().height()})};

	const double padding = 16;

	// include bounds given by vertices
	for(const Vertex* vertex : m_elems_vertices)
	{
		QPointF vertexpos = vertex->scenePos();

		if(vertexpos.x() < pt1.x())
			pt1.setX(vertexpos.x() -  padding);
		if(vertexpos.x() > pt2.x())
			pt2.setX(vertexpos.x() +  padding);
		if(vertexpos.y() < pt1.y())
			pt1.setY(vertexpos.y() -  padding);
		if(vertexpos.y() > pt2.y())
			pt2.setY(vertexpos.y() +  padding);
	}

	setSceneRect(QRectF{pt1, pt2});
}



void VisView::AddVertex(const QPointF& pos)
{
	Vertex *vertex = new Vertex{pos};
	m_elems_vertices.push_back(vertex);
	m_scene->addItem(vertex);
}


void VisView::mousePressEvent(QMouseEvent *evt)
{
	QPoint posVP = evt->pos();
	QPointF posScene = mapToScene(posVP);

	QList<QGraphicsItem*> items = this->items(posVP);
	QGraphicsItem* item = nullptr;
	bool item_is_vertex = false;

	for(int itemidx=0; itemidx<items.size(); ++itemidx)
	{
		item = items[itemidx];
		auto iter = std::find(m_elems_vertices.begin(), m_elems_vertices.end(), static_cast<Vertex*>(item));
		item_is_vertex = (iter != m_elems_vertices.end());
		if(item_is_vertex)
			break;
	}

	// only select vertices
	if(!item_is_vertex)
		item = nullptr;


	if(evt->button() == Qt::LeftButton)
	{
		// if no vertex is at this position, create a new one
		if(!item)
		{
			AddVertex(posScene);
			m_dragging = true;
			UpdateAll();
		}

		else
		{
			// vertex is being dragged
			if(item_is_vertex)
			{
				m_dragging = true;
			}
		}
	}
	else if(evt->button() == Qt::RightButton)
	{
		// if a vertex is at this position, remove it
		if(item && item_is_vertex)
		{
			m_scene->removeItem(item);
			auto iter = std::find(m_elems_vertices.begin(), m_elems_vertices.end(), static_cast<Vertex*>(item));
			if(iter != m_elems_vertices.end())
				iter = m_elems_vertices.erase(iter);
			delete item;
			UpdateAll();
		}
	}

	QGraphicsView::mousePressEvent(evt);
}


void VisView::mouseReleaseEvent(QMouseEvent *evt)
{
	if(evt->button() == Qt::LeftButton)
		m_dragging = false;

	UpdateAll();
	QGraphicsView::mouseReleaseEvent(evt);
}


void VisView::mouseMoveEvent(QMouseEvent *evt)
{
	QGraphicsView::mouseMoveEvent(evt);

	if(m_dragging)
	{
		QResizeEvent evt{size(), size()};
		resizeEvent(&evt);
		UpdateAll();
	}

	QPoint posVP = evt->pos();
	QPointF posScene = mapToScene(posVP);
	emit SignalMouseCoordinates(posScene.x(), posScene.y());
}


void VisView::ClearVertices()
{
	for(Vertex* vertex : m_elems_vertices)
	{
		m_scene->removeItem(vertex);
		delete vertex;
	}
	m_elems_vertices.clear();

	UpdateAll();
}


void VisView::UpdateAll()
{
	// get vertices
	m_vertices.clear();
	m_vertices.reserve(m_elems_vertices.size());
	std::transform(m_elems_vertices.begin(), m_elems_vertices.end(), std::back_inserter(m_vertices),
		[](const Vertex* vert) -> t_vec { return m::create<t_vec>({vert->x(), vert->y()}); } );

	if(m_sortvertices)
		std::tie(m_vertices, std::ignore) = sort_vertices_by_angle<t_vec>(m_vertices);

	UpdateEdges();
	UpdateKer();
}


void VisView::UpdateEdges()
{
	// remove previous edges
	for(QGraphicsItem* item : m_elems_edges)
	{
		m_scene->removeItem(item);
		delete item;
	}
	m_elems_edges.clear();


	QPen penEdge;
	penEdge.setStyle(Qt::SolidLine);
	penEdge.setWidthF(2.);
	penEdge.setColor(QColor::fromRgbF(0., 0., 1.));


	for(std::size_t vertidx = 0; vertidx < m_vertices.size(); ++vertidx)
	{
		std::size_t vertidx2 = (vertidx+1) % m_vertices.size();

		const t_vec& vertex1 = m_vertices[vertidx];
		const t_vec& vertex2 = m_vertices[vertidx2];

		QLineF line{QPointF{vertex1[0], vertex1[1]}, QPointF{vertex2[0], vertex2[1]}};
		QGraphicsItem *item = m_scene->addLine(line, penEdge);
		m_elems_edges.push_back(item);
	}
}


void VisView::UpdateKer()
{
	// remove previous vis poly
	for(QGraphicsItem* item : m_elems_ker)
	{
		m_scene->removeItem(item);
		delete item;
	}
	m_elems_ker.clear();


	std::vector<t_vec> verts_reversed;
	verts_reversed.reserve(m_vertices.size());
	for(auto iter=m_vertices.rbegin(); iter!=m_vertices.rend(); ++iter)
		verts_reversed.push_back(*iter);

	//auto kerpoly = calc_ker<t_vec>(m_vertices, g_eps);
	//auto kerpoly_reversed = calc_ker<t_vec>(verts_reversed, g_eps);
	auto kerpoly = calc_ker_ineff<t_vec>(m_vertices, g_eps);
	auto kerpoly_reversed = calc_ker_ineff<t_vec>(verts_reversed, g_eps);

	// in case the vertices were inserted in reversed order
	if(kerpoly_reversed.size() > kerpoly.size())
		kerpoly = kerpoly_reversed;


	QPen penKer;
	penKer.setStyle(Qt::SolidLine);
	penKer.setWidthF(2.);
	penKer.setColor(QColor::fromRgbF(1., 0., 0., 1.));

	QBrush brushKer;
	brushKer.setColor(QColor::fromRgbF(1., 0., 0., 0.1));
	brushKer.setStyle(Qt::SolidPattern);

	QPolygonF poly;
	for(std::size_t vertidx = 0; vertidx < kerpoly.size(); ++vertidx)
	{
		//std::size_t vertidx2 = (vertidx+1) % kerpoly.size();
		const t_vec& vertex1 = kerpoly[vertidx];
		//const t_vec& vertex2 = kerpoly[vertidx2];

		//QLineF line{QPointF{vertex1[0], vertex1[1]}, QPointF{vertex2[0], vertex2[1]}};
		//QGraphicsItem *item = m_scene->addLine(line, penKer);
		//m_elems_ker.push_back(item);

		poly << QPointF(vertex1[0], vertex1[1]);
	}

	QGraphicsItem *item = m_scene->addPolygon(poly, penKer, brushKer);
	m_elems_ker.push_back(item);
}


void VisView::SetSortVertices(bool b)
{
	m_sortvertices = b;
	UpdateAll();
}


// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------

VisWnd::VisWnd(QWidget* pParent) : QMainWindow{pParent},
	m_scene{new QGraphicsScene{this}},
	m_view{new VisView{m_scene.get(), this}},
	m_statusLabel{std::make_shared<QLabel>(this)}
{
	// ------------------------------------------------------------------------
	// restore settings
	QSettings settings{this};

	if(settings.contains("wnd_geo"))
	{
		QByteArray arr{settings.value("wnd_geo").toByteArray()};
		this->restoreGeometry(arr);
	}
	else
	{
		resize(1024, 768);
	}
	if(settings.contains("wnd_state"))
	{
		QByteArray arr{settings.value("wnd_state").toByteArray()};
		this->restoreState(arr);
	}

	m_view->SetSortVertices(
		settings.value("sort_vertices", m_view->GetSortVertices()).toBool());
	// ------------------------------------------------------------------------


	m_view->setRenderHints(QPainter::Antialiasing);

	setWindowTitle("Visibility");
	setCentralWidget(m_view.get());

	QStatusBar *statusBar = new QStatusBar{this};
	statusBar->addPermanentWidget(m_statusLabel.get(), 1);
	setStatusBar(statusBar);


	// menu actions
	QAction *actionNew = new QAction{"New", this};
	connect(actionNew, &QAction::triggered, [this]() { m_view->ClearVertices(); });

	QAction *actionLoad = new QAction{"Load...", this};
	connect(actionLoad, &QAction::triggered, [this]()
	{
		if(QString file = QFileDialog::getOpenFileName(this, "Load Data", "",
			"XML Files (*.xml);;All Files (* *.*)"); file!="")
		{
			std::ifstream ifstr(file.toStdString());
			if(!ifstr)
			{
				QMessageBox::critical(this, "Error", "File could not be opened for loading.");
				return;
			}

			m_view->ClearVertices();

			ptree::ptree prop{};
			ptree::read_xml(ifstr, prop);

			std::size_t vertidx = 0;
			while(true)
			{
				std::ostringstream ostrVert;
				ostrVert << "vis2d.vertices." << vertidx;

				auto vertprop = prop.get_child_optional(ostrVert.str());
				if(!vertprop)
					break;

				auto vertx = vertprop->get_optional<t_real>("<xmlattr>.x");
				auto verty = vertprop->get_optional<t_real>("<xmlattr>.y");

				if(!vertx || !verty)
					break;

				m_view->AddVertex(QPointF{*vertx, *verty});

				++vertidx;
			}

			if(vertidx > 0)
				m_view->UpdateAll();
			else
				QMessageBox::warning(this, "Warning", "File contains no data.");
		}
	});

	QAction *actionSaveAs = new QAction{"Save as...", this};
	connect(actionSaveAs, &QAction::triggered, [this]()
	{
		if(QString file = QFileDialog::getSaveFileName(this, "Save Data", "",
			"XML Files (*.xml);;All Files (* *.*)"); file!="")
		{
			std::ofstream ofstr(file.toStdString());
			if(!ofstr)
			{
				QMessageBox::critical(this, "Error", "File could not be opened for saving.");
				return;
			}

			ptree::ptree prop{};

			std::size_t vertidx = 0;
			for(const Vertex* vertex : m_view->GetVertexElems())
			{
				QPointF vertexpos = vertex->scenePos();

				std::ostringstream ostrX, ostrY;
				ostrX << "vis2d.vertices." << vertidx << ".<xmlattr>.x";
				ostrY << "vis2d.vertices." << vertidx << ".<xmlattr>.y";

				prop.put<t_real>(ostrX.str(), vertexpos.x());
				prop.put<t_real>(ostrY.str(), vertexpos.y());

				++vertidx;
			}

			ptree::write_xml(ofstr, prop, ptree::xml_writer_make_settings('\t', 1, std::string{"utf-8"}));
		}
	});

	QAction *actionExportSvg = new QAction{"Export SVG...", this};
	connect(actionExportSvg, &QAction::triggered, [this]()
	{
		if(QString file = QFileDialog::getSaveFileName(this, "Export SVG", "",
			"SVG Files (*.svg);;All Files (* *.*)"); file!="")
		{
			QSvgGenerator svggen;
			svggen.setSize(QSize{width(), height()});
			svggen.setFileName(file);

			QPainter paint(&svggen);
			m_scene->render(&paint);
		}
	});

	QAction *actionQuit = new QAction{"Exit", this};
	connect(actionQuit, &QAction::triggered, [this]()
		{ this->close(); });


	QAction *actionSort = new QAction{"Sort Vertices", this};
	actionSort->setCheckable(true);
	actionSort->setChecked(m_view->GetSortVertices());
	connect(actionSort, &QAction::toggled, [this](bool b) { m_view->SetSortVertices(b); });


	// menu
	QMenu *menuFile = new QMenu{"File", this};
	QMenu *menuCalc = new QMenu{"Calculate", this};

	menuFile->addAction(actionNew);
	menuFile->addSeparator();
	menuFile->addAction(actionLoad);
	menuFile->addAction(actionSaveAs);
	menuFile->addSeparator();
	menuFile->addAction(actionExportSvg);
	menuFile->addSeparator();
	menuFile->addAction(actionQuit);

	menuCalc->addAction(actionSort);


	// menu bar
	QMenuBar *menuBar = new QMenuBar{this};
	menuBar->setNativeMenuBar(false);
	menuBar->addMenu(menuFile);
	menuBar->addMenu(menuCalc);
	setMenuBar(menuBar);


	// connections
	connect(m_view.get(), &VisView::SignalMouseCoordinates, [this](double x, double y) -> void
	{
		SetStatusMessage(QString("x=%1, y=%2.").arg(x, 5).arg(y, 5));
	});


	SetStatusMessage("Ready.");
}


void VisWnd::SetStatusMessage(const QString& msg)
{
	m_statusLabel->setText(msg);
}


void VisWnd::closeEvent(QCloseEvent *e)
{
	// ------------------------------------------------------------------------
	// save settings
	QSettings settings{this};

	QByteArray geo{this->saveGeometry()}, state{this->saveState()};
	settings.setValue("wnd_geo", geo);
	settings.setValue("wnd_state", state);
	settings.setValue("sort_vertices", m_view->GetSortVertices());
	// ------------------------------------------------------------------------

	QMainWindow::closeEvent(e);
}


VisWnd::~VisWnd()
{
}


// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------

static inline void set_locales()
{
	std::ios_base::sync_with_stdio(false);

	::setlocale(LC_ALL, "C");
	std::locale::global(std::locale("C"));
	QLocale::setDefault(QLocale::C);
}


int main(int argc, char** argv)
{
	try
	{
		auto app = std::make_unique<QApplication>(argc, argv);
		app->setOrganizationName("tw");
		app->setApplicationName("vis2d");
		set_locales();

		auto vis = std::make_unique<VisWnd>();
		vis->show();

		return app->exec();
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}

	return -1;
}
// ----------------------------------------------------------------------------
