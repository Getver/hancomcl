const rowsPerPage = 20;
const pageButtonsToShow = 5;

const API_BASE = '/api';
let serverSort = { key: 'CHR', dir: 'asc' };

let allData = [];
let fileHeaders = [];
let currentPage = parseInt(localStorage.getItem('currentPage')) || 1;
let columnDescriptions = {};
let filters = { gene: '', snp: '', distance: '', pValueMin: '', pValueMax: '', fdrMin: '', fdrMax: '' };

const tableHead = document.getElementById('data-table-head');
const tableBody = document.getElementById('data-table-body');
const paginationControls = document.getElementById('pagination-controls');

// --------------------- 컬럼 설명 ---------------------
async function loadColumnDescriptions() {
  try {
    const res = await fetch(`${API_BASE}/column-info`);
    columnDescriptions = await res.json();
  } catch (err) {
    console.error("컬럼 설명 불러오기 실패:", err);
  }
}

// --------------------- 서버에서 페이지 가져오기 ---------------------
async function fetchEqtlPage(page = 1) {
  const params = new URLSearchParams();
  params.set('page', page);
  params.set('size', rowsPerPage);
  params.set('sort', serverSort.key);
  params.set('order', serverSort.dir);

  if (filters.gene) params.set('gene', filters.gene);
  if (filters.snp) params.set('snp', filters.snp);
  if (filters.distance) params.set('distance', filters.distance);
  if (filters.pValueMin) params.set('pValueMin', filters.pValueMin);
  if (filters.pValueMax) params.set('pValueMax', filters.pValueMax);
  if (filters.fdrMin) params.set('fdrMin', filters.fdrMin);
  if (filters.fdrMax) params.set('fdrMax', filters.fdrMax);

  const quick = document.getElementById('quick-search')?.value.trim();
  if (quick) params.set('search', quick);

  const res = await fetch(`${API_BASE}/eqtl?${params.toString()}`);
  if (!res.ok) {
    console.error("데이터 조회 실패");
    return;
  }
  const pageData = await res.json();

  fileHeaders = pageData.headers;
  allData = pageData.items;

  renderTableHead();
  renderTableRowsFromServer(pageData.total, page);
}

// --------------------- 테이블 헤더 ---------------------
function renderTableHead() {
  tableHead.innerHTML = '';
  const tr = document.createElement('tr');

  fileHeaders.forEach(h => {
    const th = document.createElement('th');
    if (columnDescriptions[h]) th.title = columnDescriptions[h]; // 툴팁
    th.textContent = h;
    th.addEventListener('click', () => sortColumn(h));
    tr.appendChild(th);
  });

  const viewTh = document.createElement('th');
  viewTh.textContent = 'VIEW';
  tr.appendChild(viewTh);
  tableHead.appendChild(tr);
}

// --------------------- 테이블 본문 ---------------------
function renderTableRowsFromServer(totalRows, page) {
  currentPage = page;
  localStorage.setItem('currentPage', currentPage);

  tableBody.innerHTML = '';
  const frag = document.createDocumentFragment();

  allData.forEach(row => {
    const tr = document.createElement('tr');
    fileHeaders.forEach(h => {
      const td = document.createElement('td');
      td.textContent = row[h] ?? '-';
      tr.appendChild(td);
    });
    tr.appendChild(createViewImages(row));  // VIEW 아이콘
    frag.appendChild(tr);
  });

  tableBody.appendChild(frag);
  updatePaginationUI_Server(totalRows);
}

// --------------------- 페이지네이션 ---------------------
function updatePaginationUI_Server(totalRows) {
  paginationControls.innerHTML = '';
  const pageCount = Math.ceil(totalRows / rowsPerPage);
  if (pageCount <= 1) return;

  const currentGroup = Math.floor((currentPage - 1) / pageButtonsToShow);
  const startPage = currentGroup * pageButtonsToShow + 1;
  const endPage = Math.min(startPage + pageButtonsToShow - 1, pageCount);

  paginationControls.appendChild(createBtnServer('«', 1, currentPage === 1));
  paginationControls.appendChild(createBtnServer('‹', currentPage - 1, currentPage === 1));
  for (let i = startPage; i <= endPage; i++) {
    const btn = createBtnServer(i, i, false);
    if (i === currentPage) btn.classList.add('active');
    paginationControls.appendChild(btn);
  }
  paginationControls.appendChild(createBtnServer('›', currentPage + 1, currentPage === pageCount));
  paginationControls.appendChild(createBtnServer('»', pageCount, currentPage === pageCount));
}

function createBtnServer(text, page, disabled) {
  const btn = document.createElement('button');
  btn.textContent = text;
  btn.disabled = disabled;
  btn.addEventListener('click', () => fetchEqtlPage(page));
  return btn;
}

// --------------------- 정렬 ---------------------
function sortColumn(header) {
  const asc = !(serverSort.key === header && serverSort.dir === 'asc');
  serverSort = { key: header, dir: asc ? 'asc' : 'desc' };

  document.querySelectorAll('.results-table thead th').forEach(th => {
    th.classList.remove('sorted-asc', 'sorted-desc');
    if (th.textContent.trim() === header) th.classList.add(asc ? 'sorted-asc' : 'sorted-desc');
  });

  fetchEqtlPage(1);
}

// --------------------- VIEW 아이콘 ---------------------
function createViewImages(row) {
  const imgInfos = [
    { src: '/image/red_arrow.png',   key: 'GENE', urlPrefix: 'https://www.gtexportal.org/home/gene/' },
    { src: '/image/green_arrow.png', key: 'GENE', urlPrefix: 'https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=' },
    { src: '/image/blue_arrow.png',  key: 'SNP',  urlPrefix: 'https://www.ncbi.nlm.nih.gov/snp/' },
    { src: '/image/orange_arrow.png',key: 'SNP',  urlPrefix: 'https://www.ebi.ac.uk/gwas/variants/' }
  ];
  const td = document.createElement('td');
  imgInfos.forEach(info => {
    const img = document.createElement('img');
    img.src = info.src;
    img.style.cssText = 'width:24px;height:24px;margin-right:4px;cursor:pointer;';
    img.addEventListener('click', () => {
      let value = row[info.key];
      if (!value) return;
      window.open(info.urlPrefix + value, '_blank');
    });
    td.appendChild(img);
  });
  return td;
}

// --------------------- CSV 다운로드 ---------------------
document.querySelector('.csv-button').addEventListener('click', async () => {
  const params = new URLSearchParams({ page: 1, size: 100000, sort: serverSort.key, order: serverSort.dir });
  if (filters.gene) params.set('gene', filters.gene);
  if (filters.snp) params.set('snp', filters.snp);
  if (filters.distance) params.set('distance', filters.distance);
  if (filters.pValueMin) params.set('pValueMin', filters.pValueMin);
  if (filters.pValueMax) params.set('pValueMax', filters.pValueMax);
  if (filters.fdrMin) params.set('fdrMin', filters.fdrMin);
  if (filters.fdrMax) params.set('fdrMax', filters.fdrMax);
  const quickVal = document.getElementById('quick-search')?.value.trim();
  if (quickVal) params.set('search', quickVal);

  const res = await fetch(`${API_BASE}/eqtl?${params.toString()}`);
  if (!res.ok) {
    alert('CSV 생성용 데이터 조회 실패');
    return;
  }
  const data = await res.json();
  if (!data.items?.length) {
    alert('다운로드할 데이터가 없습니다.');
    return;
  }

  const headers = [...data.headers, 'VIEW'];
  const csvRows = [headers.join(',')];
  data.items.forEach(row => {
    csvRows.push([...data.headers.map(h => JSON.stringify(row[h] ?? '')), ''].join(','));
  });

  const blob = new Blob([csvRows.join('\n')], { type: 'text/csv;charset=utf-8;' });
  const link = document.createElement('a');
  link.href = URL.createObjectURL(blob);
  link.download = 'snp_data.csv';
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
});

// --------------------- 이벤트 ---------------------
const quickEl = document.getElementById('quick-search');
if (quickEl) quickEl.addEventListener('input', () => fetchEqtlPage(1));

document.querySelector('.search-action-button').addEventListener('click', e => {
  e.preventDefault();
  filters.gene = document.getElementById('gene_id').value.trim();
  filters.snp = document.getElementById('snp_id').value.trim();
  filters.distance = document.getElementById('distance').value.trim();

  const pMin = document.getElementById('pval_min').value.trim();
  const pMax = document.getElementById('pval_max').value.trim();
  filters.pValueMin = pMin || '';
  filters.pValueMax = pMax || '';

  const fdrMin = document.getElementById('fdr_min').value.trim();
  const fdrMax = document.getElementById('fdr_max').value.trim();
  filters.fdrMin = fdrMin || '';
  filters.fdrMax = fdrMax || '';

  fetchEqtlPage(1);
});

document.getElementById('reset-filters').addEventListener('click', () => {
  document.querySelectorAll('.filters-panel form input').forEach(input => (input.value = ''));
  if (quickEl) quickEl.value = '';
  filters = { gene: '', snp: '', distance: '', pValueMin: '', pValueMax: '', fdrMin: '', fdrMax: '' };
  fetchEqtlPage(1);
});

// --------------------- 초기화 ---------------------
document.addEventListener('DOMContentLoaded', async () => {
  await loadColumnDescriptions();
  fetchEqtlPage(currentPage);
});
